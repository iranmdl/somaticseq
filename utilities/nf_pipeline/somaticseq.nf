#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.sample = ["tumor": ["name": "sample0.tumor", "bam": "../test-data/samples/sample0/tumor/sample0.tumor.cleaned.bam"], "normal": ["name": "sample0.normal", "bam": "../test-data/samples/sample0/normal/sample0.normal.cleaned.bam"] ]
params.regions = "../test-data/samples/sample0/regions.bed"
params.reference = "../test-data/small-resources/ref_renamed.fa"
params.indel_known = "../test-data/small-resources/known_renamed.vcf.gz"
params.dbsnp = "../test-data/small-resources/dbsnp_renamed.vcf.gz"
params.cosmic = "../test-data/small-resources/cosmic_renamed.vcf.gz"
params.classifier = ["snv": "/isilon/Analysis/datainsights/projects/TCGA_training_pairs_4LUAD_5LUSC/M2JSDULPK/sSNV_9WES.tsv.ntChange.Classifier.RData", "indel": "../test-data/small-resources/sINDEL_9WES.tsv.ntChange.Classifier.RData"]
params.truth = ["snv": "", "indel": ""]

bwaIndex = "${params.reference}.{amb,ann,bwt,pac,sa}"

params.containers = ["gatk": "broadinstitute/gatk3:3.7-0", "strelka": "lethalfang/strelka:2.7.1", "samtools": "comics/samtools", "vardict": "lethalfang/vardictjava:1.5.1", "somaticsniper": "lethalfang/somaticsniper:1.0.5.0", "scalpel": "lethalfang/scalpel:0.5.3", "mutect2": "broadinstitute/gatk:4.beta.3", "muse": "marghoob/muse:1.0rc_c", "lofreq": "marghoob/lofreq:2.1.2", "jointsnvmix2": "lethalfang/jointsnvmix2:0.7.5", "somaticseq": "lethalfang/somaticseq:2.3.0", "bwa": "biocontainers/bwa"]
params.gatk = "/usr/GenomeAnalysisTK.jar"

ref = params.reference
refDict = ref.take(ref.lastIndexOf('.')) + ".dict"

params.align = ["enabled": true]

params.markduplicates = ["enabled": false]
params.realigner = ["enabled": false]
params.baserecalibration = ["enabled": false]

params.strelka = ["enabled": true]
params.vardict = ["enabled": true, "vaf": 0.1]
params.somaticsniper = ["enabled": true, "mq": 25, "bq": 15, "prior": 0.0001]
params.scalpel = ["enabled": true]
params.mutect2 = ["enabled": true]
params.muse = ["enabled": true]
params.lofreq = ["enabled": true]
params.jointsnvmix2 = ["enabled": true, "convergence_threshold": 0.01, "skip_size": 200 ]

params.somaticseq = ["enabled": true]

dummyFileChannel = Channel.value(file("none"))
dummyPairFileChannel = Channel.value([file("none"), file("none")])

if (params.sample.tumor.bam && params.align.enabled) {
    extractFastq = true
    alignReads = true
} else {
    extractFastq = false
    alignReads = false
}

if (params.sample.tumor.fastqs) {
    extractFastq = false
    alignReads = true
}

if (extractFastq) {
    bamsForFastq = Channel.from(["normal", params.sample.normal.name, file(params.sample.normal.bam)], ["tumor", params.sample.tumor.name, file(params.sample.tumor.bam)])

    process sortBamReadname {
        tag { sampleName }
        container params.containers.samtools

        input:
        set val(sampleType), val(sampleName), file(inBam) from bamsForFastq

        output:
        set val(sampleType), val(sampleName), file('namesorted.bam') into sortedBamsByName

        """
        samtools sort -n -T sorted -o namesorted.bam $inBam
        """
    }

    process bamToFastq {
        tag { sampleName }
        container params.containers.samtools

        input:
        set val(sampleType), val(sampleName), file('namesorted.bam') from sortedBamsByName

        output:
        set val(sampleType), val(sampleName), file('reads1.fq'), file('reads2.fq') into readsForAlignment

        """
        samtools fastq -1 reads1.fq -2 reads2.fq namesorted.bam
        """
    }
} else if (alignReads) {
    readsForAlignment = Channel.from(["normal", params.sample.normal.name, file(params.sample.normal.fastqs[0]), file(params.sample.normal.fastqs[1])], ["tumor", params.sample.tumor.name, file(params.sample.tumor.fastqs[0]), file(params.sample.tumor.fastqs[1])])
} 

 
if (alignReads) {
    process bwaAlign {
        tag { sampleName }
        validExitStatus 0
        container params.containers.bwa

        input:
        set val(sampleType), val(sampleName), file(reads1), file(reads2) from readsForAlignment
        set file(referenceFile), file("*") from Channel.value([file(params.reference), file(bwaIndex)])

        output:
        set val(sampleType), val(sampleName), file('aligned.sam') into alignedSams

        script:
        """
        #!/usr/bin/env bash

        set -o pipefail -e
        bwa mem -p -R \"@RG\tID:${sampleName}\tSM:${sampleName}\tLB:${sampleName}\tPL:Illumina\" $referenceFile ${reads1} ${reads2} | grep -v "^@PG" > aligned.sam
        """
    }

    process sortSamCoordinateIndex {
        tag { sampleName }
        container params.containers.samtools

        input:
        set val(sampleType), val(sampleName), file(alignedSam) from alignedSams

        output:
        set val(sampleType), val(sampleName), file("${sampleType}.sorted.bam"), file("${sampleType}.sorted.bam.bai") into sortedBamsByCoordinate

        """
        samtools view -1 -F 2304 -o aligned.bam ${alignedSam}
        samtools sort -T sorted -o ${sampleType}.sorted.bam aligned.bam && samtools index ${sampleType}.sorted.bam
        """
    }
} else {
    sortedBamsByCoordinate = Channel.from(["normal", params.sample.normal.name, file(params.sample.normal.bam), file(params.sample.normal.bam + ".bai")], ["tumor", params.sample.tumor.name, file(params.sample.tumor.bam), file(params.sample.tumor.bam + ".bai")])
}

if (params.markduplicates.enabled) {
    process markDuplicates {
        tag { sampleName }
        container params.containers.samtools

        input:
        set val(sampleType), val(sampleName), file("in.bam"), file("in.bam.bai") from sortedBamsByCoordinate

        output:
        set val(sampleType), val(sampleName), file("${sampleType}.dupmarked.bam"), file("${sampleType}.dupmarked.bam.bai") into dupmarkedBams

        script:
        """
        samtools rmdup in.bam ${sampleType}.dupmarked.bam && samtools index ${sampleType}.dupmarked.bam
        """
    }
} else {
    sortedBamsByCoordinate.set { dupmarkedBams }
}


if (params.realigner.enabled) {
    dupmarkedBams.reduce( [:] ) { a, b -> a[b[0]] = b.subList(1, 4); return a }
                 .into { bamsForTargetCreator; bamsForRealign }

    process gatkRealignerTargetCreator {
        container params.containers.gatk

        input:
        val bams from bamsForTargetCreator
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict), file(indelKnown), file(indelKnownIndex) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict), file(params.indel_known), file(params.indel_known + ".tbi")])

        output:
        file "indel.intervals" into realignIntervals

        script:
        tumorBam = bams["tumor"][1]
        normalBam = bams["normal"][1]

        """
        java -jar ${params.gatk} -T RealignerTargetCreator -R $referenceFile --known $indelKnown -I ${tumorBam} -I ${normalBam} -o indel.intervals
        """
    }

    process gatkRealigner {
        container params.containers.gatk

        input:
        val bams from bamsForRealign
        file realignIntervals
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])

        output:
        set val("tumor"), val(params.sample.tumor.name), file("tumor.dupmarked.JointRealigned.bam") into tumorRealigned
        set val("normal"), val(params.sample.normal.name), file("normal.dupmarked.JointRealigned.bam") into normalRealigned

        script:
        tumorBam = bams["tumor"][1]
        normalBam = bams["normal"][1]

        """
        java -jar ${params.gatk} -T IndelRealigner -nWayOut .JointRealigned.bam -R $referenceFile -targetIntervals ${realignIntervals} -I ${tumorBam} -I ${normalBam}
        """
    }

    tumorRealigned.mix(normalRealigned).set { realignedBams }
} else {
    dupmarkedBams.set { realignedBams } 
}



if (params.baserecalibration.enabled) {

    realignedBams.into { bamsForBaseRecalibrator; bamsForPrintReads }

    process gatkBaseRecalibrator {
        tag { sampleName }
        container params.containers.gatk

        input:
        set val(sampleType), val(sampleName), file("in.bam") from bamsForBaseRecalibrator
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict), file(dbsnp), file(dbsnpIndex) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict), file(params.dbsnp), file(params.dbsnp + ".tbi")])

        output:
        file "${sampleType}.recal.table" into recalTables

        script:
        """
        java -jar ${params.gatk} -T BaseRecalibrator -R $referenceFile -knownSites $dbsnp -I in.bam -o ${sampleType}.recal.table
        """
    }

    process gatkPrintReads {
        tag { sampleName }
        container params.containers.gatk
    
        input:
        set val(sampleType), val(sampleName), file("in.bam") from bamsForPrintReads
        file recalTables
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])

        output:
        set val(sampleType), val(sampleName), file("${sampleType}.recalibrated.bam"), file("${sampleType}.recalibrated.bam.bai") into recalibratedBams

        script:
        """
        java -jar ${params.gatk} -T PrintReads -R $referenceFile -I in.bam -BQSR ${recalTables} -o ${sampleType}.recalibrated.bam
        """
    }

    recalibratedBams.set { bamsForSomaticCalling }
} else {
    realignedBams.set { bamsForSomaticCalling }
}

process indexBams {
    tag { sampleName }
    container params.containers.samtools
    
    input:
    set val(sampleType), val(sampleName), file(inBam), file('dummy.bai') from bamsForSomaticCalling

    output:
    set val(sampleType), file("${inBam.getName()}"), file("${inBam.getName()}.bai") into indexedBams

    script:
    """
    samtools index $inBam
    """
}
// Now run the multiple somatic calling tools


indexedBams.reduce( [:] ) { a, b -> a[b[0]] = b.subList(1, b.size()); return a }
           .map { [it.tumor[0], it.tumor[1], it.normal[0], it.normal[1]] }
           .into { bamsForStrelka; bamsForMutect; bamsForScalpel; bamsForJointSNVMix; bamsForLofreq; bamsForVarDict; bamsForMuse; bamsForSomaticSniper; bamsForSomaticSeq }

if (params.strelka.enabled) {
    process runStrelka {
        container params.containers.strelka

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForStrelka
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])

        output:
        set file("results/variants/somatic.snvs.vcf.gz"), file("results/variants/somatic.indels.vcf.gz") into strelkaOut

        script:
        """
        /opt/strelka2/bin/configureStrelkaSomaticWorkflow.py --tumorBam $tumorBam --normalBam $normalBam --referenceFasta $referenceFile --callMemMb=4096 --runDir .
        ./runWorkflow.py -m local -j 1
        """
    }
} else {
    strelkaOut = dummyPairFileChannel
}

if (params.vardict.enabled) {
    process runVarDict {
        container params.containers.vardict

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForStrelka
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
        file regions from file(params.regions)

        output:
        file "vardict.vcf" into vardictOut

        script: 
        """
        /opt/VarDict-1.5.1/bin/VarDict -G $referenceFile -b \"${tumorBam}|${normalBam}\" -Q 1 -c 1 -S 2 -E 3 -g 4 -f ${params.vardict.vaf} -h ${regions}> vardict.var
        cat vardict.var | awk 'NR!=1' | /opt/VarDict/testsomatic.R | /opt/VarDict/var2vcf_paired.pl -N 'TUMOR|NORMAL' -f ${params.vardict.vaf} > vardict.vcf
        """
    }
} else {
    vardictOut = dummyFileChannel
}

if (params.somaticsniper.enabled) {
    process runSomaticSniper {
        container params.containers.somaticsniper

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForSomaticSniper
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
 
        output:
        file "somaticsniper.vcf" into somaticsniperOut

        script:
        """
        /opt/somatic-sniper/build/bin/bam-somaticsniper -q ${params.somaticsniper.mq} -Q ${params.somaticsniper.bq} -s ${params.somaticsniper.prior} -F vcf -f ${referenceFile} ${tumorBam} ${normalBam} somaticsniper.vcf
        """
    }
} else {
    somaticsniperOut = dummyFileChannel
}

if (params.scalpel.enabled) {
    process runScalpel {
        container params.containers.scalpel

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForScalpel
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
        file regions from file(params.regions)

        output:
        file "scalpel.vcf" into scalpelOut

        script:
        """
        /opt/scalpel-0.5.3/scalpel-discovery --somatic --ref ${referenceFile} --bed ${regions} --normal ${normalBam} --tumor ${tumorBam} --window 600 --dir .
        /opt/scalpel-0.5.3/scalpel-export --somatic --db main/somatic.db.dir --ref ${referenceFile} --bed ${regions} > scalpel_tmp.vcf
        cat scalpel_tmp.vcf | /opt/vcfsorter.pl ${referenceFileDict} - > scalpel.vcf
        """
    }
} else {
    scalpelOut = dummyFileChannel
}

if (params.mutect2.enabled) {
    process runMutect2 {
        container params.containers.mutect2

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForMutect
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
        file regions from file(params.regions)
        set file(dbsnp), file(dbsnpIndex) from Channel.value([file(params.dbsnp), file(params.dbsnp + ".tbi")])

        output:
        file "mutect2.vcf" into mutectOut

        script:
        """
        java -Xmx8g -jar /gatk/build/libs/gatk.jar Mutect2 --reference ${referenceFile} --intervals ${regions} --dbsnp ${dbsnp} --input ${tumorBam} --input ${normalBam} --normalSampleName ${params.sample.normal.name} --tumorSampleName ${params.sample.tumor.name} --output mutect_unfiltered.vcf
        java -Xmx8g -jar /gatk/build/libs/gatk.jar FilterMutectCalls --variant mutect_unfiltered.vcf --output mutect2.vcf
        """
    }
} else {
    mutectOut = dummyFileChannel
}

if (params.muse.enabled) {
    process runMuse {
        container params.containers.muse

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForMuse
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
        file regions from file(params.regions)
        set file(dbsnp), file(dbsnpIndex) from Channel.value([file(params.dbsnp), file(params.dbsnp + ".tbi")])

        output:
        file "muse.vcf" into museOut

        script:
        """
        awk '{print \$1 \"\t\" \$2 \"\t\" \$3}' ${regions} > bed_3columns.bed
        MuSEv1.0rc_submission_c039ffa call -O muse -l ${regions} -f ${referenceFile} ${tumorBam} ${normalBam}
        MuSEv1.0rc_submission_c039ffa sump -I muse.MuSE.txt -E -O muse.vcf -D ${dbsnp}
        """
    }
} else {
    museOut = dummyFileChannel
}

if (params.lofreq.enabled) {
    process runLofreq {
        container params.containers.lofreq

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForLofreq
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
        file regions from file(params.regions)
        set file(dbsnp), file(dbsnpIndex) from Channel.value([file(params.dbsnp), file(params.dbsnp + ".tbi")])

        output:
        set file("lofreqsomatic_final.snvs.vcf.gz"), file("lofreqsomatic_final.indels.vcf.gz") into lofreqOut

        script:
        """
        lofreq somatic -t ${tumorBam} -n ${normalBam} --call-indels -l ${regions} -f ${referenceFile} -d ${dbsnp} -o lofreq
        """
    }
} else {
    lofreqOut = dummyPairFileChannel
}

if (params.jointsnvmix2.enabled) {
    process runJointSNVMix2 {
        container params.containers.jointsnvmix2

        input:
        set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForJointSNVMix
        set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
        file regions from file(params.regions)
        set file(dbsnp), file(dbsnpIndex) from Channel.value([file(params.dbsnp), file(params.dbsnp + ".tbi")])

        output:
        file "JointSNVMix2.vcf" into jointsnvmix2Out
        
        script:
        """
        /opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py train joint_snv_mix_two --convergence_threshold ${params.jointsnvmix2.convergence_threshold} --skip_size ${params.jointsnvmix2.skip_size} ${referenceFile} ${normalBam} ${tumorBam} /opt/JointSNVMix-0.7.5/config/joint_priors.cfg /opt/JointSNVMix-0.7.5/config/joint_params.cfg jsm.parameter.cfg

        echo -e '##fileformat=VCFv4.1' > JointSNVMix2.vcf
        echo -e '##INFO=<ID=AAAB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and AB in Tumor\">' >> JointSNVMix2.vcf
        echo -e '##INFO=<ID=AABB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and BB in Tumor\">' >> JointSNVMix2.vcf
        echo -e '##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">' >> JointSNVMix2.vcf
        echo -e '##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">' >> JointSNVMix2.vcf
        echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR' >> JointSNVMix2.vcf

        /opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py classify joint_snv_mix_two ${referenceFile} ${normalBam} ${tumorBam} jsm.parameter.cfg /dev/stdout | awk 'NR!=1 && \$4!=\"N\" && \$10+\$11>=0.95' | awk '{print \$1 \"\t\" \$2 \"\t.\t\" \$3 \"\t\" \$4 \"\t.\t.\tAAAB=\" \$10 \";AABB=\" \$11 \"\tRD:AD\t\" \$5 \":\" \$6 \"\t\" \$7 \":\" \$8}' | /opt/vcfsorter.pl ${referenceFileDict} - >> JointSNVMix2.vcf
        """
    }
} else {
    jointsnvmix2Out = dummyFileChannel
}

process runSomaticSeq {
    container params.containers.somaticseq

    input:
    file jointsnvmix2Variants from jointsnvmix2Out
    file mutect2Variants from mutectOut
    file somaticsniperVariants from somaticsniperOut
    file vardictVariants from vardictOut
    file museVariants from museOut
    file scalpelVariants from scalpelOut
    set file(lofreqSNVs), file(lofreqIndels) from lofreqOut
    set file(strelkaSNVs), file(strelkaIndels) from strelkaOut
    set file(dbsnp), file(dbsnpIndex) from Channel.value([file(params.dbsnp), file(params.dbsnp + ".tbi")])
    set file(cosmic), file(cosmicIndex) from Channel.value([file(params.cosmic), file(params.cosmic + ".tbi")])
    file snvClassifier from Channel.value(file(params.classifier.snv))
    file indelClassifier from Channel.value(file(params.classifier.indel))
    set file(referenceFile), file(referenceFileIndex), file(referenceFileDict) from Channel.value([file(params.reference), file(params.reference + ".fai"), file(refDict)])
    set file(tumorBam), file(tumorBai), file(normalBam), file(normalBai) from bamsForSomaticSeq

    when:
    params.somaticseq.enabled

    script:
    mutect2Text = params.mutect2.enabled ? "--mutect2 ${mutect2Variants}" : ""
    jsmText = params.jointsnvmix2.enabled ? "--jsm ${jointsnvmix2Variants}" : ""
    vardictText = params.vardict.enabled ? "--vardict ${vardictVariants}" : ""
    sniperText = params.somaticsniper.enabled ? "--sniper ${somaticsniperVariants}" : ""
    museText = params.muse.enabled ? "--muse ${museVariants}" : ""
    scalpelText = params.scalpel.enabled ? "--scalpel ${scalpelVariants}": ""
    strelkaSNVText = params.strelka.enabled ? "--strelka-snv ${strelkaSNVs}": ""
    strelkaIndelText = params.strelka.enabled ? "--strelka-indel ${strelkaIndels}": ""
    lofreqSNVText = params.lofreq.enabled ? "--lofreq-snv ${lofreqSNVs}": ""
    lofreqIndelText = params.lofreq.enabled ? "--lofreq-indel ${lofreqIndels}": ""
    
    """
    /opt/somaticseq/SomaticSeq.Wrapper.sh --output-dir out --genome-reference ${referenceFile} --tumor-bam ${tumorBam} --normal-bam ${normalBam} ${jsmText} ${mutect2Text} ${vardictText} ${sniperText} ${museText} ${scalpelText} ${strelkaSNVText} ${strelkaIndelText} ${lofreqSNVText} ${lofreqIndelText} --dbsnp ${dbsnp} --cosmic ${cosmic} --gatk /opt/GATK/GenomeAnalysisTK.jar --ada-r-script /opt/somaticseq/r_scripts/ada_model_predictor.R --classifier-snv ${snvClassifier} --classifier-indel ${indelClassifier}
    """
}


workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.cawMessage()
  log.info "Workflow execution stopped with the following message: " + workflow.errorMessage
}
