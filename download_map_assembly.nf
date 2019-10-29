#!/usr/bin/env nextflow

params.assemblies_gz = "*_genomic.fna.gz"
params.projects_accs = "Project_accessions_head.txt"
params.publish = "demo"
params.alienseq = "$baseDir/databases//alienTrimmerPF8contaminants.fasta"
params.minlength = 45
params.cpus=4
params.mismatch=1
compressed_assemblies = file(params.assemblies_gz)
project_accessions = file(params.projects_accs)
myDir = file(params.publish)
myDir.mkdirs()


/*
 * Step 1. Builds the genome index required by the mapping process. Concatenate, and gunzip first.
 */
process buildIndex {
    publishDir "$myDir/index", mode: 'copy'

    input:
    file assemblies_gz from compressed_assemblies
      
    output:
    file("assemblies.index*") into assemblies_index
        
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
	//publishDir "$myDir/demo_${x}_runAccs" 
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
	publishDir "$myDir/demo_${x}_reads", mode: 'copy'

	input:
	val x from project2
	val run
	
	output:
    // if no pair we don't take this project
	set val(x), file('*_1.fastq'), file('*_2.fastq')  optional true into readChannel
    //val '*.fastq' into reads

    shell:
	"""
	fasterq-dump -e 2 !{run}
	"""
}



/*process qc {
	publishDir "$myDir/demo_${pair_id}_qc"
	
	input:
	set pair_id, file(forward), file(reverse) from readChannel

	output:
	file "fastqc_${pair_id}_logs" into fastqc_ch

	shell:
	"""
	mkdir fastqc_!{pair_id}_logs
	fastqc -o fastqc_!{pair_id}_logs !{forward} !{reverse} --quiet
	"""

}*/
/*Channel.from(reads)
    .map { it -> tuple(it.split('_')[0], file(it))}
    .groupTuple()
    .set{ readChannel }*/

//Channel.from(readChannel)
//       .subscribe { println it }
       //.map { file -> tuple(file.baseName, file) }
       //.groupTuple()
       
//.splitText(sep: "_")

process filtering {
    //publishDir "$myDir", mode: 'copy'
    cpus params.cpus
    memory "8G"

    input:
    set pair_id, file(forward), file(reverse) from readChannel
    file(leptospira_genome) from assemblies_index

    output:
    set pair_id, file("mapped/*.1"), file("mapped/*.2") into mappedChannel

    shell:
    """
    mkdir mapped
    bowtie2 -q -N !{params.mismatch} -1 !{forward} -2 !{reverse} \
            -x assemblies.index --al-conc mapped/ -S /dev/null \
            -p !{params.cpus} --very-sensitive-local
    """
}


process trim {
	tag "${pair_id}" 

	input:
	set pair_id, file(forward), file(reverse) from mappedChannel

	output:
	set pair_id, file("*_R1.fastq"), file("*_R2.fastq") into trimChannel

	shell:
	"""
	AlienTrimmer -if ${forward} -ir ${reverse} -of ${pair_id}_R1.fastq \
                 -or ${pair_id}_R2.fastq -os ${pair_id}_sgl.fastq \
                 -c ${params.alienseq} -l ${params.minlength}
	"""
}

process assembly {
    publishDir "$myDir", mode: 'copy'
    memory "30G"
    cpus params.cpus
    tag "${pair_id}" 

    input:
    //set pair_id, file(forward), file(reverse) from khmerChannel
    set pair_id, file(forward), file(reverse) from trimChannel

    output:
    set pair_id, file("assembly/*_{spades}.fasta") into contigsChannel

    shell:
    """
    mkdir assembly
    spades.py -1 !{forward} -2 !{reverse} -t !{params.cpus} -o assembly/
    rename_fasta.py -i assembly/scaffolds.fasta \
            -s !{pair_id} -o assembly/!{pair_id}_spades.fasta
    """
}
