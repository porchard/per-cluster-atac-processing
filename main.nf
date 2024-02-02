#!/usr/bin/env nextflow

GENOME = params.genome
MARKERS = params.markers

nextflow.enable.dsl=2

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
	'hg38': (1..22).collect({it -> 'chr' + it}),
	'rn5': (1..20).collect({it -> 'chr' + it}),
	'rn6': (1..20).collect({it -> 'chr' + it}),
	'mm9': (1..19).collect({it -> 'chr' + it}),
	'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
	'hg38': 'human',
	'rn5': 'rat',
	'rn6': 'rat',
	'mm9': 'mouse',
	'mm10': 'mouse'
]

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]

def get_macs2_genome_size (genome) {
	return MACS2_GENOME_SIZE[genome]
}

def get_blacklists (genome) {
    if (params.blacklist[genome] instanceof String) {
        return [params.blacklist[genome]]
    } else {
        return params.blacklist[genome]
    }
}

def make_excluded_regions_arg (genome) {
    return get_blacklists(genome).collect({'--excluded-region-file ' + it}).join(' ')
}

def get_genome_size  (genome) {
	MACS2_GENOME_SIZE[genome]
}

def get_tss (genome) {
	params.tss[genome]
}

def get_organism (genome) {
	ORGANISMS[genome]
}

def get_chrom_sizes (genome) {
	params.chrom_sizes[genome]
}


process aggregate_bam {

	publishDir "${params.results}/bam/per-library-pass-qc-nuclei"
	container "library://porchard/default/general:20220107"

	input:
	tuple val(library), path("to-subset.bam"), path(clusters)

	output:
	tuple val(library), path("${library}.bam"), path("${library}.bam.bai")

	"""
	convert-barcode.py --rna-barcodes ${params.rna_barcodes} --atac-barcodes ${params.atac_barcodes} --from RNA --to ATAC --field 1 $clusters | cut -f1 > keep-nuclei.atac.txt
	filter-bam-by-barcode.py to-subset.bam ${library}.bam keep-nuclei.atac.txt
	samtools index ${library}.bam
	"""

}


process per_library_per_cluster_bam {

	publishDir "${params.results}/bam/per-library-per-cluster"
	container "library://porchard/default/general:20220107"
	maxForks 10

	input:
	tuple val(library), path("bam-to-split"), path(clusters)

	output:
	tuple val(library), path("*.bam")
	
	"""
	convert-barcode.py --rna-barcodes ${params.rna_barcodes} --atac-barcodes ${params.atac_barcodes} --from RNA --to ATAC --field 1 $clusters > barcode_to_cluster.atac.txt
    split-bam-by-cluster.py bam-to-split ${library}. barcode_to_cluster.atac.txt
	"""

}


process merge_cluster {

	publishDir "${params.results}/bam/per-cluster"
	container "library://porchard/default/general:20220107"
	maxForks 5

	input:
	tuple val(cluster), path(bams)

	output:
	tuple val(cluster), path("${cluster}.bam")

	"""
	samtools merge ${cluster}.bam ${bams.join(' ')}
	"""

}


process make_aggregate {

	publishDir "${params.results}/bam/aggregate"
	container "library://porchard/default/general:20220107"
	maxForks 5

	input:
	path(bams)

	output:
	tuple val('all'), path("all.bam")

	"""
	samtools merge all.bam ${bams.join(' ')}
	"""

}


process bamtobed {

	tag "${cluster}"
	container "library://porchard/default/general:20220107"
	maxForks 5

	input:
	tuple val(cluster), path(bam)

	output:
	tuple val(cluster), path("reads.bed")

	"""
	bedtools bamtobed -i $bam > reads.bed
	"""

}


process call_peaks {

	publishDir "${params.results}/peaks/broad"
	container "library://porchard/default/general:20220107"
	tag "${cluster}"
	memory '10 GB'

	input:
	tuple val(cluster), path(reads)

	output:
	tuple val(cluster), path("${cluster}_treat_pileup.bdg"), emit: bedgraph
	tuple val(cluster), path("${cluster}_peaks.broadPeak.noblacklist"), emit: peaks

	"""
	macs2 callpeak -t $reads --outdir . -f BED -n ${cluster} --SPMR -g ${get_genome_size(GENOME)} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
	bedtools intersect -a ${cluster}_peaks.broadPeak -b ${get_blacklists(GENOME).join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}_peaks.broadPeak.noblacklist
	"""

}


process call_narrow_peaks {

	publishDir "${params.results}/peaks/narrow"
	container "library://porchard/default/general:20220107"
	tag "${cluster}"
	memory '10 GB'

	input:
	tuple val(cluster), path(reads)

	output:
	tuple val(cluster), path("${cluster}_peaks*"), path("${cluster}_summits*")
	tuple val(cluster), path("${cluster}_summits.bed"), emit: summits

	"""
    macs2 callpeak -t $reads --qvalue 0.05 --outdir . -f BED -n ${cluster} --SPMR -g ${get_genome_size(GENOME)} --nomodel --shift -37 --seed 762873 --extsize 73 -B --keep-dup all --call-summits
	bedtools intersect -a ${cluster}_summits.bed -b ${get_blacklists(GENOME).join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}_summits.noblacklist
	bedtools intersect -a ${cluster}_peaks.narrowPeak -b ${get_blacklists(GENOME).join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}_peaks.narrowPeak.noblacklist
	"""

}


process lafave_summit_processing {

	publishDir "${params.results}/peaks/summit-extension"
	container "library://porchard/default/general:20220107"
	memory '10 GB'

	input:
	tuple val(cluster), path(summits)

	output:
	tuple val(cluster), path("${cluster}_peaks.noblacklist")

	"""
	cat $summits | awk '\$5>=3' > summits-fdr-filtered.bed
	lafave-summit-processing.py --extension 150 summits-fdr-filtered.bed > ${cluster}_peaks.bed
	bedtools intersect -a ${cluster}_peaks.bed -b ${get_blacklists(GENOME).join(' ')} -v | sort -k1,1 -k2n,2 > ${cluster}_peaks.noblacklist
	"""

}


process cluster_bigwigs {

	publishDir "${params.results}/bigwig"
	container "library://porchard/default/general:20220107"
	memory '50 GB'
	tag "${cluster}"

	input:
	tuple val(cluster), path(bedgraph)

	output:
	tuple val(cluster), path("${cluster}.bw")

	"""
	bedClip $bedgraph ${get_chrom_sizes(GENOME)} bedgraph.clipped.bdg
	LC_COLLATE=C sort -k1,1 -k2n,2 bedgraph.clipped.bdg > bedgraph.sorted.bdg
	bedGraphToBigWig bedgraph.sorted.bdg ${get_chrom_sizes(GENOME)} ${cluster}.bw
	"""
	
}


process plot_marker_genes {

	publishDir "${params.results}/plot-marker-gene-signal"
	container "library://porchard/default/general:20220107"
	memory '20 GB'

	input:
	path(bigwigs)
	val(markers)
	path(gene_bed)

	output:
	path("per-cluster-atac.pdf")

	"""
	plot-atac-signal-at-marker-genes.py --gene-bed $gene_bed --genes ${markers.split(',').join(' ')} --bigwigs ${bigwigs.join(' ')} --out per-cluster-atac.pdf
	"""

}


workflow {

	libraries = params.libraries.keySet()
	
	input_bams = []
	clusters = []
    
	for (library in libraries) {
		input_bams << [library, file(params.libraries[library].bam)]
		clusters << [library, file(params.libraries[library].clusters)]
    }

	input_bams = Channel.from(input_bams) // library, bam
	clusters = Channel.from(clusters) // library, clusters
	
	input_bams.combine(clusters, by: 0) | aggregate_bam
	bams_per_library_per_cluster = (input_bams.combine(clusters, by: 0) | per_library_per_cluster_bam).transpose().map({it -> [it[0], it[1].getName().tokenize('.')[1], it[1]]}) // library, cluster, bam
	bams_per_cluster = bams_per_library_per_cluster.map({it -> [it[1], it[2]]}).groupTuple() | merge_cluster // cluster, bam
	bam_aggregate = bams_per_cluster.map({it -> it[1]}).toSortedList() | make_aggregate // 'all', bam
	read_beds = bams_per_cluster.mix(bam_aggregate) | bamtobed
	peaks = call_peaks(read_beds) // cluster, peaks
	summits = call_narrow_peaks(read_beds).summits | lafave_summit_processing
	
	bigwigs = peaks.bedgraph | cluster_bigwigs
	plot_marker_genes(bigwigs.map({it -> it[1]}).toSortedList(), MARKERS, Channel.fromPath(params.gene_bed[GENOME]))
}
