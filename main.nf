#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
Metagenomic Taxonomy Pipeline:
    
Takes an input of fastq files (zipped or unzipped), performs a de novo assembly and identifies species via blastx and blastn,
producing tally files for these runs, extracted fasta files and a heatmap of viruses or other microorganisms in the host sample

USAGE: nextflow run main.nf [options] --input INPUT_DIR --output OUTPUT_DIR 

OPTIONS:

--input INPUT_DIR - [Required] A directory containing paired-end fastq files

--output OUTPUT_DIR - [Required] A directory to place output files (If not existing, pipeline will create)

--local_blast_nt BLAST_DATABASE_DR - [Required] A directory containing the blast nt database

--local_diamond DIAMOND_DIR - [Required] A directory containing a local diamond database

OPTIONAL:

    --host_fasta HOST_FASTA - A host fasta file that the reads will be aligned to to remove host contamination.
            
    --host_bt2_index INDEX_DIRECTORY - A directory containing an existing bowtie2 index to be used for host read removal. Must be in its own directory. If using a large genome, this option will greatly improve the pipeline runtime.

    --threads INT - The number of threads that can be use to run pipeline tools in parallel

    --ref REFERENCE_FASTA - The pipeline will align contigs produced by assembly to this reference

    --minLen INT - The minimum length of a read to keep post trimming [Default = 75bp]

    --minTrimQual INT - The average basecall quality threshold below which to trim a read. During trimming, trimmomatic performs a sliding window checking the average base quality, and removing the rest of the read if it drops below this treshold. [Default = 20]
    
    --ncbi_tax_dir - Directory containing ncbi taxonomy database to set up Taxonomizr 

    --blast_tax_dir - Directory containing ncbi taxonomy database locally to make blast taxonomically aware
    """
}

// Function that checks whether a directory ends in a trailing slash.
// This is useful for directory variables that are not parsed into
// file objects in the pipeline (such as the output directory).
def checkDirectoryEnding (fileName) {
    // Grabs the last character in the directory name.
    lastChar = fileName.substring(fileName.length() - 1, fileName.length())

    // Checks whether the last character is slash
    if (lastChar != "/") {
        // If it is not a slash, add that to the directory name.
        fileName = fileName + "/"
    }
    
    // Return the directory name.
    return fileName
}

// Function creates the header for the summary file based on the parameters
// supplied by the user. Because the pipeline dynamically changes based on
// what the user specifies, the summary file must also be alter to reflect
// the analysis. This will also make incorporating new modules easier.
def createSummaryHeader (hostRef, hostIdx) {
    
    // The header will always start with the sample.
    FinalHeader = 'Sample,'

    // The summary sheet will also always contain the Raw and trimmed read counts.
    FinalHeader = FinalHeader + "Raw Reads,Trimmed Reads,Deduped Reads,"

    // Next, the user may supply host removal, which would be the next useful
    // statistic to know. Thus, if a host reference or bowtie2 index is supplied,
    // add this field to the header.
    if (hostRef != false || hostIdx != false) {
        FinalHeader = FinalHeader + "Non-Host Reads,"
    }


    // Finally, the pipeline will always report the number of contigs and scaffolds
    FinalHeader = FinalHeader + "Contigs Generated,Scaffolds Generated"

    return FinalHeader
}

// If the help parameter is supplied, link display the help message
// and quit the pipeline
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// Defines input parameters. Setting to false by default
// allows us to check that these have been set by the user.
params.input = false
params.host_fasta = false
params.host_bt2_index = false
params.ref = false
params.output = false
params.always_trim_3p_bases = false
params.always_trim_5p_bases = false
params.minLen = 75
params.minTrimQual = 20
params.single_read = false
params.alignmentMode = "--local"
params.phred = 33
params.classify_singletons = false
params.max_blast_nt_evalue = "1e-10"
params.max_blastx_nr_evalue = "1e-3"
params.tally_filter_kingdom = false
//memory parameters
params.threads = 10
params.scripts_bindir = false
params.minimum_contig_length = 200
//params for databases
params.blast_tax_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
params.local_blast_nt = false
params.local_diamond = false
params.blast_tax_dir = false
params.ncbi_tax_dir = false
params.ncbi_tax_db = "ncbi_taxonomy.sqlite3"




// Inports modules
include { Setup } from "./modules.nf"
include { Index_Host_Reference } from './modules.nf'
include { QC_Report } from './modules.nf'
// The same module cannot be used more than once,
// thus it is aliased to be used multiple times.
include { QC_Report as QC_Report_Trimmed } from './modules.nf'
include { QC_Report as QC_Report_Deduped } from './modules.nf'
include { QC_Report as QC_Report_Host_Removed } from './modules.nf'
include { Trimming } from './modules.nf'
include { Remove_PCR_Duplicates } from './modules.nf'
include { Host_Read_Removal } from './modules.nf'
include { Spades_Assembly } from './modules.nf'
include { Retrieve_Contigs} from './modules.nf'
include { Quantify_Read_Mapping} from './modules.nf'
//include { Merge_Contigs_and_Singletons} from './modules.nf'
include { Contig_Alignment } from './modules.nf'
include { Write_Summary } from './modules.nf'
include { Ncbi_Tax_Setup} from './modules.nf'
include { Check_Blast_Tax} from './modules.nf'
include { Blastn_Contigs} from './modules.nf'
include { Process_Blastn_Output} from './modules.nf'
include { Tally_Blastn_Results} from './modules.nf'
include { Distribute_Blastn_Results} from './modules.nf'
include { Blastx_Remaining_Contigs} from './modules.nf'
include { Split_Merged_Blastx_Results} from './modules.nf'
include { Tally_Blastx_Results} from './modules.nf'
include { Distribute_Blastx_Results} from './modules.nf'
include { Virus_Mapping_Matrix } from './modules.nf'
include { Create_Heatmap } from './modules.nf'
adapters = file("${baseDir}/adapters.fa")

// Checks the input parameter
inDir = ""
if (params.input == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No input directory provided. Pipeline requires an input directory."
    exit(1)
}
else if (!(file(params.input).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.input} is not an existing directory."
    exit(1)
}
else {
    inDir = file(params.input).toString()
}

println "Input Directory: ${inDir}"

// Create a channel for the input files.

if (params.single_read != false) {
    // Single-end reads
    inputFiles_ch = Channel.fromPath("${inDir}/*.f*q*")
        .map { file -> 
            baseName = file.getBaseName() // Extract the base name of the file
            [baseName, file] // Create a tuple with the base name and the file
        }
} else {
    // Paired-end reads
    inputFiles_ch = Channel.fromFilePairs("${inDir}/*.f*q*") { file ->
        // Extract everything before the first underscore
        file.getBaseName().split("_")[0]
    }
}
// Checks the output parameter.
outDir = ''
if (params.output == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No output directory provided. Pipeline requires an output directory."
    exit(1)
}
else {
    // If the parameter is set, ensure that the directory provided ends
    // in a trailing slash (to keep things consistent throughout) the
    // pipeline code.
    outDir = file(params.output).toString()
}
println outDir

// Parses the host options (--host_fasta and --host_bt2_index).
// For this, we cannot use a channel like we did for the input files
// because then Nextfow would only run other modules once.
// Thus, we need to manually create a tuple of input data to pass to the indexing
// step and the alignment step.
hostRefData = ''
hostRefIdxData = ''
hostRefName = 'NONE'
if (params.host_fasta != false && params.host_bt2_index != false) {
    // If both options are supplied, notify the user and exit.
    println "ERROR: you have specified both a host fasta file and bowtie2 index. Please only supply one."
    exit(1)
}
else if (params.host_fasta != false) {
    if (!(file(params.host_fasta).exists())) {
        // If the file supplied does not exist, notify the user and exit.
        println "ERROR: ${params.host_fasta} does not exist."
        exit(1)
    }
    else {
        // Parse the file into a file object
        hostRef = file(params.host_fasta)
        // Use the getBaseName() function to 
        // get the reference name. This will be
        // used to name the bowtie2 index.
        hostRefName = hostRef.getBaseName()
        // Place these both into a tuple.
        hostRefData = tuple(hostRefName, hostRef)
    }
}
// If the user supplied the --host_bt2_index
else if (params.host_bt2_index != false) {
    if (!(file(params.host_bt2_index).exists())) {
        // If the index provided does not exist, notify the user and exit.
        println "Error: ${params.host_bt2_index} does not exist."
        exit(1)
    }
    else {
        // Parse the directory into a file object
        hostRefDir = file(params.host_bt2_index)
        println hostRefDir
        // Grab a list of file objects from the directory
        // ending in .bt2
        indexFiles = file("${hostRefDir}/*.bt2")
        if (indexFiles.size() == 0){
            // If there are no file in the directory ending in bt2, notify the user and exit
            println "Index Directory provided (${params.host_bt2_index}) does not contain any bt2 files"
            exit(1)
        }
        else {
            // Use the getSimpleName() function to grab the base name
            // of the index files (getBaseName() removes anything following
            // the first . in a file name.)
            hostRefName = indexFiles[0].getSimpleName()
            println hostRefName
            // Place the index dir and name into a tuple.
            hostRefIdxData = tuple(hostRefDir, hostRefName)
        }
    }
}


// Parses the ref option.
refFile = ''
if (params.ref != false) {
    if (!(file(params.ref).exists())) {
        // If the reference file did not exist, notify the user and exit.
        println "ERROR: ${params.ref} does not exist."
        exit(1)
    }
    else {
        // Parse the provided file into a file object.
        refFile = file(params.ref)
    }
}

//checking for blast databases
blast_nt_Dir = ""
if (params.local_blast_nt == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No blast database directory provided. Pipeline requires an blast database directory."
    exit(1)
}
else if (!(file(params.local_blast_nt).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.local_blast_nt} is not an existing directory."
    exit(1)
}
else {
    blast_nt_Dir = file(params.local_blast_nt).toString()
}

//create channel for local nt database
blast_nt_database = Channel.value("${blast_nt_Dir}/nt")


//checks for local diamond
diamond_Dir = ""
if (params.local_diamond == false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No diamond database directory provided. Pipeline requires an diamond database directory."
    exit(1)
}
else if (!(file(params.local_diamond).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.local_diamond} is not an existing directory."
    exit(1)
}
else {
    diamond_Dir = file(params.local_diamond).toString()
}
diamond_database = Channel.value("${diamond_Dir}/nr.dmnd")
    
//check for local blast tax directory path
blast_tax_ch = Channel.empty()

// if this path was provided as a parameter, then create a channel
// from this path and set a boolean to true to indicate it's an existing
// directory
if (params.blast_tax_dir) {
   blast_tax_ch = Channel.fromPath( params.blast_tax_dir )
                         .map { path -> [ path , true ] }  
} else {
   // if this path was *not* provided as a parameter, then create a channel
   // from a bogus path and set a boolean to false 
   // to indicate it *doesn't* refer to an existing directory
   blast_tax_ch = Channel.fromPath( "does_not_exist" )
                         .map { path -> [ path , false ] }  
}
//check for filter option on tally step
//TODO: currently not working; maybe try and fix
/*
if (params.tally_filter_kingdom) {
   filter = Channel.of( params.tally_filter_kingdom )
                         .map { val -> [ val , true ] }  
} else {
   // if this path was *not* provided as a parameter, then create a channel
   // from a bogus path and set a boolean to false 
   // to indicate it *doesn't* refer to an existing directory
   filter = Channel.of( "does_not_exist" )
                         .map { val -> [ val , false ] }  
}
*/

//check for ncbi directory output path
ncbi_tax_Dir = ""
if (params.ncbi_tax_dir== false) {
    // If the parameter is not set, notify the user and exit.
    println "ERROR: No blast database directory provided. Pipeline requires an blast database directory."
    exit(1)
}
else if (!(file(params.ncbi_tax_dir).isDirectory())) {
    // If the input directory is not set, notify the user and exit.
    println "ERROR: ${params.ncbi_tax_dir} is not an existing directory."
    exit(1)
}
else {
    ncbi_tax_Dir = file(params.ncbi_tax_dir).toString()
}



println "Input Directory: ${inDir}"

assembler = 'SPADES'

summaryHeader = createSummaryHeader(params.host_fasta, params.host_bt2_index)

total_deduped = 0

workflow {

    Setup( summaryHeader, hostRefName, params.minLen, params.minTrimQual, outDir )

    // Use FASTQC to perform an initial QC check on the reads
    QC_Report( inputFiles_ch, outDir, "FASTQC-Pre-Processing", params.threads )

    // Perform adapter and quality trimming with Trimmomatic.
    Trimming( inputFiles_ch, outDir, adapters, params.minLen, params.minTrimQual)

    // Use FASTQC to perform a QC check on the trimmed reads.
    QC_Report_Trimmed( Trimming.out[0], outDir, "FASTQC-Trimmed", params.threads )

    // Perform PCR Duplicate removal using prinseq.
    Remove_PCR_Duplicates( Trimming.out[0], outDir, Trimming.out[2] , total_deduped)

    // Use FASTQC to perform a QC check on the deduped reads.
    QC_Report_Deduped( Remove_PCR_Duplicates.out[0], outDir, "FASTQC-Deduplicated", params.threads )

    //uses taxonomizr to preprocess and create sqlite database for ncbi taxonomy
    Ncbi_Tax_Setup(ncbi_tax_Dir)
    setup_ncbi_dir = Ncbi_Tax_Setup.out[0].collect()

    //check blast for local taxdb and if taxonomically aware
    Check_Blast_Tax(blast_tax_ch)
    checked_blast_tax = Check_Blast_Tax.out[0].collect()
    

    // If the user supplied a reference fasta file, build a bowtie2 index and use
    // that for alignment.
    if (params.host_fasta) {
        Index_Host_Reference( hostRefData, outDir, params.threads )
        Host_Read_Removal( Remove_PCR_Duplicates.out[0], outDir, Index_Host_Reference.out, params.alignmentMode, params.threads, Remove_PCR_Duplicates.out[1] )
        QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQC-Host-Removed", params.threads )

        
        // Perform de novo assembly using spades.
        Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, params.phred, Host_Read_Removal.out[2] )
        if (params.ref != false) {
            // Align the contigs to a reference genome using minimap2 and samtools
            Contig_Alignment( Spades_Assembly.out[0], outDir, refFile )
            
            
        }

        Write_Summary( Spades_Assembly.out[2], outDir )
        //retrieve contigs from spades to conver to 1 line
        Retrieve_Contigs( Spades_Assembly.out[0], Spades_Assembly.out[3], outDir)
        //creates contig weights file
        Quantify_Read_Mapping(Retrieve_Contigs.out[0], Retrieve_Contigs.out [1], params.threads, outDir)
        //merge contigs & singletons
        //Merge_Contigs_and_Singletons(Quantify_Read_Mapping.out[0], Quantify_Read_Mapping.out[2])
        //use blast nt database sesarch to first find matches on nucleotide level
        Blastn_Contigs(Spades_Assembly.out[0], checked_blast_tax, params.max_blast_nt_evalue, blast_nt_database, params.threads, outDir)
        //adds header to blast data & outputs fasta file of contigs/singletons that didn't match
        Process_Blastn_Output(Blastn_Contigs.out[0], Quantify_Read_Mapping.out[1], outDir)
        //perform tally on blast results
        Tally_Blastn_Results(Process_Blastn_Output.out[2], setup_ncbi_dir, outDir, Remove_PCR_Duplicates.out[2])
        //create seperate fasta files for top hits of blastn search
        Distribute_Blastn_Results(Process_Blastn_Output.out[1], setup_ncbi_dir, outDir, Tally_Blastn_Results.out[1])

        //collect all unmatched files after blastn search into one file to enter for diamond
        Process_Blastn_Output.out[3]
                .collectFile(name:'merged_contigs.fasta')
                .set{merged_blastx_ch}
                
        //use blastx/diamond to search nr database for contigs that didn't match previously
           
        Blastx_Remaining_Contigs(merged_blastx_ch, diamond_database, checked_blast_tax, params.threads)

        //combine blastx output with rest of info from blastn
        Blastx_Remaining_Contigs.out[0]
            .combine(Process_Blastn_Output.out[0])     
            .set{merged_blastx_with_input_ch}
        //split back into respective files    
        Split_Merged_Blastx_Results(merged_blastx_with_input_ch, outDir)
        //tally results from blastx
        Tally_Blastx_Results(Split_Merged_Blastx_Results.out[0], setup_ncbi_dir, outDir, Remove_PCR_Duplicates.out[2])
        //distribute into fasta files from results
        Distribute_Blastx_Results(Split_Merged_Blastx_Results.out[1], setup_ncbi_dir, outDir, Tally_Blastn_Results.out[1])
        //combine all tally files to one channel
        Tally_Blastx_Results.out[0]
            .mix(Tally_Blastn_Results.out[0])
            .collect()
            .set{all_tally_ch}
        //create mapping matrix for taxas
        Virus_Mapping_Matrix(all_tally_ch, outDir)
        //generate heatmap from mapping matrix
        Create_Heatmap(Virus_Mapping_Matrix.out[0], outDir)
    
    }
    // If the user supplied an existing bowtie2 index, use that for alignment.
    else if (params.host_bt2_index) {
        Host_Read_Removal( Remove_PCR_Duplicates.out[0], outDir, hostRefIdxData, params.alignmentMode, params.threads, Remove_PCR_Duplicates.out[1] )
        QC_Report_Host_Removed( Host_Read_Removal.out[0], outDir, "FASTQC-Host-Removed", params.threads )
        // Perform de novo assembly using spades.
        Spades_Assembly( Host_Read_Removal.out[0], outDir, params.threads, params.phred, Host_Read_Removal.out[2] )
            if (params.ref != false) {
                // Align the contigs to a reference genome using minimap2 and samtools
                Contig_Alignment( Spades_Assembly.out[0], outDir, refFile )
            }

            Write_Summary( Spades_Assembly.out[2], outDir )
            //retrieve contigs from spades to conver to 1 line
            Retrieve_Contigs( Spades_Assembly.out[0], Spades_Assembly.out[3], outDir)
            //creates contig weights file
            Quantify_Read_Mapping(Retrieve_Contigs.out[0], Retrieve_Contigs.out [1], params.threads, outDir)
            //merge contigs & singletons
            //Merge_Contigs_and_Singletons(Quantify_Read_Mapping.out[0], Quantify_Read_Mapping.out[2])
            //use blast nt database sesarch to first find matches on nucleotide level
            Blastn_Contigs(Spades_Assembly.out[0], checked_blast_tax, params.max_blast_nt_evalue, blast_nt_database, params.threads, outDir)
            //adds header to blast data & outputs fasta file of contigs/singletons that didn't match
            Process_Blastn_Output(Blastn_Contigs.out[0], Quantify_Read_Mapping.out[1],outDir)
            //perform tally on blast results
            Tally_Blastn_Results(Process_Blastn_Output.out[2], setup_ncbi_dir, outDir, Remove_PCR_Duplicates.out[2])
            //create seperate fasta files for top hits of blastn search
            Distribute_Blastn_Results(Process_Blastn_Output.out[1], setup_ncbi_dir, outDir, Tally_Blastn_Results.out[1])

            //collect all unmatched files after blastn search into one file to enter for diamond
            Process_Blastn_Output.out[3]
                    .collectFile(name:'merged_contigs.fasta')
                    .set{merged_blastx_ch}
                    
            //use blastx/diamond to search nr database for contigs that didn't match previously
            
            Blastx_Remaining_Contigs(merged_blastx_ch, diamond_database, checked_blast_tax, params.threads)

            //combine blastx output with rest of info from blastn
            Blastx_Remaining_Contigs.out[0]
                .combine(Process_Blastn_Output.out[0])     
                .set{merged_blastx_with_input_ch}
            //split back into respective files    
            Split_Merged_Blastx_Results(merged_blastx_with_input_ch, outDir)
            //tally results from blastx
            Tally_Blastx_Results(Split_Merged_Blastx_Results.out[0], setup_ncbi_dir, outDir, Remove_PCR_Duplicates.out[2])
            //distribute into fasta files from results
            Distribute_Blastx_Results(Split_Merged_Blastx_Results.out[1], setup_ncbi_dir, outDir, Tally_Blastn_Results.out[1])
            //combine all tally files to one channel
            Tally_Blastx_Results.out[0]
                .mix(Tally_Blastn_Results.out[0])
                .collect()
                .set{all_tally_ch}
            //create mapping matrix for taxas
            Virus_Mapping_Matrix(all_tally_ch, outDir)
            //generate heatmap from mapping matrix
            Create_Heatmap(Virus_Mapping_Matrix.out[0], outDir)
    
    }
    else {    
        // Perform de novo assembly using spades.
        Spades_Assembly( Remove_PCR_Duplicates.out[0], outDir, params.threads, params.phred, Remove_PCR_Duplicates.out[1] )
        if (params.ref != false) {
            // Align the contigs to a reference genome using minimap2 and samtools
            Contig_Alignment( Spades_Assembly.out[0], outDir, refFile )
            
        }

        Write_Summary( Spades_Assembly.out[2], outDir )
        //retrieve contigs from spades to conver to 1 line
        Retrieve_Contigs( Spades_Assembly.out[0], Spades_Assembly.out[3], outDir)
        //creates contig weights file
        Quantify_Read_Mapping(Retrieve_Contigs.out[0], Retrieve_Contigs.out [1], params.threads, outDir)
        //merge contigs & singletons
        //Merge_Contigs_and_Singletons(Quantify_Read_Mapping.out[0], Quantify_Read_Mapping.out[2])
        //use blast nt database sesarch to first find matches on nucleotide level
        Blastn_Contigs(Spades_Assembly.out[0], checked_blast_tax, params.max_blast_nt_evalue, blast_nt_database, params.threads, outDir)
        //adds header to blast data & outputs fasta file of contigs/singletons that didn't match
        Process_Blastn_Output(Blastn_Contigs.out[0], Quantify_Read_Mapping.out[1],outDir)
        //perform tally on blast results
        Tally_Blastn_Results(Process_Blastn_Output.out[2], setup_ncbi_dir, outDir, Remove_PCR_Duplicates.out[2])
        //create seperate fasta files for top hits of blastn search
        Distribute_Blastn_Results(Process_Blastn_Output.out[1], setup_ncbi_dir, outDir, Tally_Blastn_Results.out[1])

        //collect all unmatched files after blastn search into one file to enter for diamond
        Process_Blastn_Output.out[3]
                .collectFile(name:'merged_contigs.fasta')
                .set{merged_blastx_ch}
                
        //use blastx/diamond to search nr database for contigs that didn't match previously
           
        Blastx_Remaining_Contigs(merged_blastx_ch, diamond_database, checked_blast_tax, params.threads)

        //combine blastx output with rest of info from blastn
        Blastx_Remaining_Contigs.out[0]
            .combine(Process_Blastn_Output.out[0])     
            .set{merged_blastx_with_input_ch}
        //split back into respective files    
        Split_Merged_Blastx_Results(merged_blastx_with_input_ch, outDir)
        //tally results from blastx
        Tally_Blastx_Results(Split_Merged_Blastx_Results.out[0], setup_ncbi_dir, outDir, Remove_PCR_Duplicates.out[2])
        //distribute into fasta files from results
        Distribute_Blastx_Results(Split_Merged_Blastx_Results.out[1], setup_ncbi_dir, outDir, Tally_Blastn_Results.out[1])
        //combine all tally files to one channel
        Tally_Blastx_Results.out[0]
            .mix(Tally_Blastn_Results.out[0])
            .collect()
            .set{all_tally_ch}
        //create mapping matrix for taxas
        Virus_Mapping_Matrix(all_tally_ch, outDir)
        //generate heatmap from mapping matrix
        Create_Heatmap(Virus_Mapping_Matrix.out[0], outDir)
    }
    
    
}