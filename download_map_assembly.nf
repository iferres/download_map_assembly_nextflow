#!/usr/bin/env nextflow

params.assemblies_gz = "*_genomic.fna.gz"
params.projects_accs = "Project_accessions_head.txt"
params.publish = "${MYSCRATCH}/demo"

compressed_assemblies = file(params.assemblies_gz)
project_accessions = file(params.projects_accs)
myDir = file(params.publish)
myDir.mkdirs()


/*
 * Step 1. Builds the genome index required by the mapping process. Concatenate, and gunzip first.
 */
process buildIndex {
    publishDir "$myDir/demo_index"

    input:
    file assemblies_gz from compressed_assemblies
      
    output:
    file 'assemblies.index*' into assemblies_index
        
    shell:
    """
    cat ${assemblies_gz} > all_assemblies_genomic.fna.gz
    gunzip all_assemblies_genomic.fna.gz
    bowtie2-build all_assemblies_genomic.fna assemblies.index
    """
}

/*
 * Takes a .txt with project accessions, and process each line
 */

Channel
	.fromPath("${params.projects_accs}")
	.splitText().map{it -> it.trim()}
	//.subscribe { println "value: $it"}
    	.set{ project }

/*
 * Downloads list of runs from project accession
 */


process obtainRunAccessions {
	publishDir "$myDir/demo_${x}_runAccs"
	
	label "onlocal"
	tag "${x}"
	
	input:
	val x from project
	
	output:
	file '*run_accessions.txt' into run_accessions
	val x into project2

	shell:
	"""
	curl 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=${x}&result=read_run&fields=run_accession' -o "${x}_run_accessions.txt"
	tail -n +2 "${x}_run_accessions.txt" > "${x}_run_accessions.txt.tmp" && mv "${x}_run_accessions.txt.tmp" "${x}_run_accessions.txt"
	"""
}


//Channel
//    .fromPath('*_run_accessions.txt')
//    .splitText(by: 1).map{it -> it.trim()}
//    .set{ runs }

run = run_accessions.splitText().map {it -> it.trim()}

process downloadReads {
	label "onlocal"
	publishDir "$myDir/demo_${x}_reads"

	input:
	val x from project2
	val run
	
	output:
	file '*_{1,2}.fastq' into reads
	val x into project3
	
	shell:
	"""
	echo !{run}
	fasterq-dump -e 1 !{run}
	"""
}


process qc {
	publishDir "$myDir/demo_${x}_qc"
	
	input:
	val x from project3
	file "" from reads

	output:
	file "fastqc_${pair_id}_logs" into fastqc_ch

	script:
	

	shell:
	"""
	mkdir fastqc_!{pair_id}_logs
	fastqc -o fastqc_!{pair_id}_logs -f fastq --quiet !{fq}
	"""

}

/*

process trim {
	publishDir 

	input:

	output:

	shell:
	"""
	"""

}


process mapReads {
  input:
  file read from reads
  file index from assemblies_index

  output:
  


}

process spadesAssembly
*/
