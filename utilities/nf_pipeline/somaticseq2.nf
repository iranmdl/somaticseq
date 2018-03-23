#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv



runextractFastq = params.runextractFastq
runfastq2bam = params.runfastq2bam
runfastqc = params.runfastqc
runMutect = params.runMutect
runMutect2 = params.runMutect2
runVarscan2 = params.runVarscan2
runSomaticsniper = params.runSomaticsniper
runVardict = params.runVardict
runMuse = params.runMuse
runJointsnvmix2 = params.runJointsnvmix2
runLofreq = params.runLofreq
runScalpel = params.runScalpel
runStrelka = params.runStrelka
runSomaticseq = params.runSomaticseq
runSomaticseq_train = params.runSomaticseq_train




// Samples file parsing
// colnames: PatID, SampleType, R1, R2, bam
// content: PatID, Normal/Tumor, abs path to fastq_r1, abs path to fastq_r2, abs path to bam if skip the alignment
sampleData = []
samples = file(params.samples)
samplesData = parseCsv(samples.text, separator: "\t")


for ( row in samplesData ) {
    row = row.toMap()
    sampleData << row
}


///////////////// extract fastq from bam ////////////////
if (runextractFastq) {
    // TODO
}

///////////////// run fastqc //////////////////
if (runfastqc){
    //TODO
}
///////////////// fastq2bam //////////////////////////////
if (runfastq2bam) {

    \\ alignment 
    process bwaAlign {
        
        tag("${PatID}_${SampleType}")

        fastqChannel = channel.from(sampleData)
        validExitStatus 0,42
        container params.containers.bwa
        cpus 24
        memory '10 GB'
        clusterOptions = '-l h_vmem=24G ' + clusterOptions

        publishDir "${finalDir}/bams/", mode: 'link'

        input:
            set val(PatID), val(SampleType), file(read1), file(read2), val(bam) from fastqChannel
            set file(ref) from params.reference

        output:
            set val(PatID), val(SampleType), file('${PatID}.${SampleType}.bam') into alignBams

        script:
            """
            bwa mem \
            -R '@RG\tID:${PatID}_${SampleType}\tLB:${PatID}_${SampleType}\tPL:illumina\tSM:${PatID}_${SampleType}' \
            -t 24 \
            ref \
            read1 \
            read2 \
            | samtools view -Sbh - \
            | samtools sort -m 8G --threads 24 -o '${PatID}.${SampleType}.bam'
            """
    }

    \\ index bam 
    process alignBamIndex {
        tag("${PatID}_${SampleType}")
        validExitStatus 0,42
        container params.containers.bwa
        cpus 1
        memory '4 GB'
        clusterOptions = '-l h_vmem=4G ' + clusterOptions

        publishDir "${finalDir}/bams/", mode: 'link'

        input:
            set val(PatID), val(SampleType), val(bam) from alignBams

        output:
            set val(PatID), val(SampleType), file('${PatID}.${SampleType}.bam') into alignBamsIndexed

        script:
            """
            samtools index bam
            """

    }

    \\ mark duplicates
    process markDup{
        tag("${PatID}_${SampleType}")
        validExitStatus 0,42
        container params.containers.picard
        cpus 1
        memory '10 GB'
        clusterOptions = '-l h_vmem=10G ' + clusterOptions

        publishDir "${finalDir}/bams/", mode: 'link'

        input:
            set val(PatID), val(SampleType), val(bam) from alignBamsIndexed

        output:
            set val(PatID), val(SampleType), file('${PatID}.${SampleType}.bam') into markdupBams

        script:
            """
            java -Xmx8g -jar /opt/picard.jar MarkDuplicatesWithMateCigar \
            I=bam \
            M=${PatID}.${SampleType}.markdup \
            CREATE_INDEX=true \
            MINIMUM_DISTANCE=1000 \
            O=${PatID}.${SampleType}.markdup.bam

            mv ${PatID}.${SampleType}.markdup.bai ${PatID}.${SampleType}.markdup.bam.bai
            """
    }

    markdupBamGroup = []
    markdupBams.inject([:]){result, e -> 
        if(! result.find{ it.key == e[0] }?.value) result[e[0]] = [:]
        result[e[0]][e[1]] = e[2]; return result}.each() { k, v -> markdupBamGroup << [k,v.Normal, v.Tumor ]} 


    \\indel realignment
    process indelRealigner{
        tag("${PatID}")
        validExitStatus 0,42
        container params.containers.gatk3
        cpus 1
        memory '10 GB'
        clusterOptions = '-l h_vmem=10G ' + clusterOptions

        publishDir "${finalDir}/bams/", mode: 'link'

        input:
            set val(PatID), file(NormalBam), file(TumorBam) from markdupBamGroup
            set file(ref) from params.reference
            set file(dbSNP) from params.dbSNP

        output:
            set val(PatID), val("Normal"), file(NormalBam.replace(".bam",".jointRealigned.bam")) into bamsForBRNormal
            set val(PatID), val("Tumor"), file(TumorBam.replace(".bam",".jointRealigned.bam")) into bamsForBRTumor

        script:
            """
            java -Xmx8g -jar GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
            -R ref \
            -I NormalBam \
            -I TumorBam \
            --known dbSNP \
            -o '${PatID}.markdup.T.N.joinRealigned.intervals'

            java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
            -T IndelRealigner \
            -R ref \
            -I NormalBam \
            -I TumorBam \
            -targetIntervals '${PatID}.markdup.T.N.joinRealigned.intervals' \
            -nWayOut .jointRealigned.bam

            mv NormalBam.replace(".bam",".jointRealigned.bai") NormalBam.replace(".bam",".jointRealigned.bam.bai")
            mv TumorBam.replace(".bam",".jointRealigned.bai") TumorBam.replace(".bam",".jointRealigned.bam.bai")
            """

    }

    bamsForBaseRecalibrator = bamsForBRNormal + bamsForBRTumor

    \\ baserecalibrator
    process indelRealigner{
        tag("${PatID}")
        validExitStatus 0,42
        container params.containers.gatk3
        cpus 1
        memory '10 GB'
        clusterOptions = '-l h_vmem=10G ' + clusterOptions

        publishDir "${finalDir}/bams/", mode: 'link'

        input:
            set val(PatID), file(SampleType), file(Bam) from bamsForBaseRecalibrator
            set file(ref) from params.reference
            set file(dbSNP) from params.dbSNP

        output:
            set val(PatID), val(SampleType), file('${PatID}.${SampleType}.markdup.jointRealigned.BQSR.bam') into finalBams

        script:
            """
            java -Xmx8g -jar GenomeAnalysisTK.jar \
            -T BaseRecalibrator \
            -R ref \
            -I Bam \
            -knownSites dbSNP \
            -o ${PatID}.${SampleType}.table

            java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
            -T PrintReads \
            -R ref \
            -I Bam \
            -BQSR ${PatID}.${SampleType}.table \
            -o '${PatID}.${SampleType}.markdup.jointRealigned.BQSR.bam'

            mv '${PatID}.${SampleType}.markdup.jointRealigned.BQSR.bai' '${PatID}.${SampleType}.markdup.jointRealigned.BQSR.bam.bai'
            """

    }    
    finalBamGroup = []
    finalBams.inject([:]){result, e -> 
        if(! result.find{ it.key == e[0] }?.value) result[e[0]] = [:]
        result[e[0]][e[1]] = e[2]; return result}.each() { k, v -> finalBamGroup << [k,v.Normal, v.Tumor ]} 
}else{
    finalBamGroup = []
    sampleData.inject([:]){result, e -> 
        if(! result.find{ it.key == e[0] }?.value) result[e[0]] = [:]
        result[e[0]][e[1]] = e[2]; return result}.each() { k, v -> finalBamGroup << [k,v.Normal, v.Tumor ]}
    
}

\\\\\\\\\\\\\\\\\\\\ Run mutation calling \\\\\\\\\\\\\\\\\\\\\\\\\\\

runMutect = params.runMutect
runMutect2 = params.runMutect2
runVarscan2 = params.runVarscan2
runSomaticsniper = params.runSomaticsniper
runVardict = params.runVardict
runMuse = params.runMuse
runJointsnvmix2 = params.runJointsnvmix2
runLofreq = params.runLofreq
runScalpel = params.runScalpel
runStrelka = params.runStrelka
runSomaticseq = params.runSomaticseq
runSomaticseq_train = params.runSomaticseq_train



\\run Mutect
if(runMutect2){
    
}


\\run Varscan2
if(runVarscan2){
    
}


\\run Somaticsniper
if(runSomaticsniper){
    
}


\\run Vardict
if(runVardict){
    
}


\\run muse
if(runMuse){
    
}


\\run Jointsnvmix2
if(runJointsnvmix2){
    
}


\\run Mutect
if(runMutect2){
    
}

