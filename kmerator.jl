#!/usr/bin/env julia

# AUTHOR : Sébastien RIQUIER, IRMB, Montpellier

print_with_color(:blue, "                                                       
                      <((((((\\\\\\
                     \/      . }\\
                     ;--..--._|}
  (\\                 '--/\--'  )
   \\\\                | '-'  :'|
    \\\\               . -==- .-|
     \\\\               \\.__.'   \\--._
     [\\\\          __.--|       //  _/'--.
     \\ \\\\       .'-._ ('-----'/ __/      \
      \\ \\\\     /   __>|      | '--.       |
       \\ \\\   |   \\   |     /    /       /
        \\ '\\ /    \\  |     |  _/       /
         \\  \\       \\ |     | /        /
          \\\\\\\  \\      \\      _  /   
 ____  __.                        | |            
|    |/ _| ___ _ __ _ __ ___  __ _| |_ ___  _ __   
|      <  | '_ ` _ \\/ _ \\ '__/ _` | __/ _ \\| '__|  
|    |  \\ | | | | ||  __/ |   (_| | || (_) | |     
|____|__ \\|_| |_| |_\\___|_|  \\__,_|\\__\\___/|_|     
        \\/                                              
-------------------------------------------------- 
Dependencies : 
- R + (stringi) & rjson libraries
- Jellyfish
- 
                                                   
                                                   
                                                  
--------------------------------------------------
")


#print_with_color(:orange,intro)


@everywhere using ParallelDataTransfer
using ArgParse
@everywhere using FastaIO
#@everywhere using ProgressMeter
using RCall

#Parse argument
s = ArgParseSettings()
@add_arg_table s begin

"--selection"
    help = "list of gene to select directly inside your fasta transcriptome file"
    action = :append_arg 
    nargs = '*'
"--appris"
    help = " indicate : 'homo_sapiens','mus_musculus','rattus_norvegicus','danio_rerio','sus_scrofa','' or virtually any specie available on APPRIS database (Rodriguez JM, Rodriguez-Rivas J, Domenico TD, Vázquez J, Valencia A, and Tress ML. Nucleic Acids Res. Database issue; 2017 Oct 23.).
    A gene have multiple possible sequence depending of its variants. this option select the principal transcript of a gene, based of APPRIS database, or by length if no data or no connection. If 'APPRIS option and not available data or no principal transcript (non-coding genes), use length instead"
    arg_type = String

"--unannotated", "-u"
      action = :store_true
      help = "activated if the provided initial fasta file correspond to an annotation external from Ensembl. Otherwise, use ensembl fasta files !"

      "--stringent", "-s"
        action = :store_true
        help = "if you think a gene-specific tag must be unique but present in ALL KNOWN TRANSCRIPT, if false, the kmer is considered as gene-specific is present only in variant of the corresponding gene, independently of the number know variant containing the kmer"

      "--verbose", "-v"
        action = :store_true
        help = "If you want this script talk too much"

    "--genome", "-g"
        help = "the genome fasta (.fa) or index by jellyfish for kmer request"
        #type = String
        required = true
    "--transcriptome", "-t"
        help = "the transcriptome fasta (.fa) (FASTA ONLY) for kmer request and transcriptional variant informations"
        #type = String
        required = true
    "--level", "-l"
        help = "type 'gene', 'transcript' or 'chimera' to extract specific kmer at these different levels. 'chimera' option must be done with 'unannotated' option!."
        required = true
    "--output", "-o"
        help = "directory of output"
        default = "."
    "--length"
        help = "length required for the kmer generation"
        default = 31
        arg_type = Int
    "fasta_file" #argument, no option
        help = "the fasta input file"
        required = true

    "--threshold"
        help = "FOR GENE LEVEL ONLY : Minimum fraction of annotated transcripts containing this kmer to admit it (default 0.5)"
        default = 0.5
        arg_type = Float64
end
parsed_args = parse_args(ARGS, s)

print(parsed_args)


# variable attribution

output = parsed_args["output"]
if ismatch(r"\/$", output)
  output = replace(output, r"/$", "")
end


println("output directory: $output")

select_option = parsed_args["selection"]
println(select_option)
if !isempty(select_option)
select_option = select_option[1]
end
APPRIS_option = parsed_args["appris"]
if APPRIS_option != nothing
    println("APPRIS selection of principal transcripts for $APPRIS_option" )
end

verbose_option = parsed_args["verbose"]
unannotated_option = parsed_args["unannotated"]



stringent_option = parsed_args["stringent"]
if stringent_option == true
println("stringent: yes")
end
admission_threshold = parsed_args["threshold"]
genome = parsed_args["genome"]
transcriptome = parsed_args["transcriptome"]
println("transcriptome : $transcriptome")
println("genome : $genome")

kmer_length = parsed_args["length"]
println("length: $kmer_length")
nbthreads =  nworkers()
println("nb process = $nbthreads")
fastafile = parsed_args["fasta_file"]
println("input sequences files = $fastafile")

level = parsed_args["level"]
println("level of kmer specificity = $level")
if unannotated_option == true
println("annotated references : no")
else
println("annotated references : yes")
end

##### Options verification part
if ! isfile("$genome")
  error("genome not found")
end
if ! isfile("$transcriptome")
  error("transcriptome not found")
end
if ! isfile("$fastafile")
  error("input fasta file not found")
end

## Ensembl annotation are ENSTXXXX.XX or ENSGXXXXX.XX in ensembl annotation, to avoid bad recognition for selection option,
## putting annotation like ENSTXXXXX.X for selection is forbidden
if select_option != nothing
    for i in select_option
        println("requested gene : $i")
        if contains(i, ".") && contains(i, "ENS")
        error("ensembl annotations with a point like ENSTXXXXXXXX.1 is forbidden, just remove the .1")
end
end
end


#println("essais")
#println(APPRIS_option != nothing)
#println(level != "gene" || unannotated_option == true)
if APPRIS_option != nothing && (level != "gene" || unannotated_option == true)
error("APPRIS option works only with the gene annotated level ")
end

#####


#verif of fasta ensembl annotation

if ismatch(r".*\.fa", transcriptome) == true
  println("input fasta (*.fa) file, continue")
else
  error("error provided transcriptome : not a fasta file") 
end

# create dictionary of transcriptome fasta
println("create dictionary of transcriptome fasta")



dico_transcriptome = Dict()

FastaReader("$transcriptome") do fr
for (desc, seq) in fr
desc_array = split(desc)
      # println( desc_array[1])
if unannotated_option == true
gene_name = desc
dico_transcriptome["$gene_name"] = seq
else
gene_name = replace(desc_array[7], "gene_symbol:", "")
ensembl_transcript_name = split(desc_array[1],'.')[1]
ensembl_gene_name = split(replace(desc_array[4], "gene:",""), '.')[1]
dico_transcriptome["$gene_name:$ensembl_transcript_name"] = seq
end

end
end

println("transcripts_distionarty finished")

#directories creation
if ! isdir("$output/sequences")
    mkdir("$output/sequences")
end
if ! isdir("$output/sequences/$kmer_length")
    mkdir("$output/sequences/$kmer_length")
end
if ! isdir("$output/tags")
    mkdir("$output/tags")
end
if ! isdir("$output/tags/$kmer_length")
    mkdir("$output/tags/$kmer_length")
end
if ! isdir("$output/jellyfish_indexes")
    mkdir("$output/jellyfish_indexes")
end
if ! isdir("$output/jellyfish_indexes/$kmer_length")
    mkdir("$output/jellyfish_indexes/$kmer_length")
end


#verif of fasta file
if ismatch(r".*\.fa", fastafile) == true
  println("input fasta (*.fa) file, continue")
else
  error("input : not a fasta file")
end

#verif of genome/transcriptome are fasta, do the jellyfish count in this case
if ismatch(r".*\.fa", genome) == true
  run(`jellyfish count -C -m $kmer_length -s 10000 -t $nbthreads -o $output/jellyfish_indexes/$kmer_length/$(replace(basename(genome), ".fa", ".jl")) $genome`)
  genome = "$output/jellyfish_indexes/$kmer_length/$(replace(basename(genome), ".fa", ".jl"))"
  println("genome kmer index output : $genome")
end
if ismatch(r".*\.fa", transcriptome) == true
  run(`jellyfish count -m $kmer_length -s 10000 -t $nbthreads -o $output/jellyfish_indexes/$kmer_length/$(replace(basename(transcriptome), ".fa", ".jl")) $transcriptome`)
  transcriptome = "$output/jellyfish_indexes/$kmer_length/$(replace(basename(transcriptome), ".fa", ".jl"))"
  println("transcriptome kmer index output : $transcriptome")
end

APPRIS_function = function(gene_ref)
#gene_ref = "ENSG00000000457"
url = "http://apprisws.bioinfo.cnio.es/rest/exporter/id/"*"$APPRIS_option"*"/" * "$gene_ref" * "?methods=appris&format=json&sc=ensembl"
println(gene_ref)
println(url)
@rput gene_ref url

R"
start_time = Sys.time()
library(stringi)
library(rjson)
#queryId <- url
DATA <- fromJSON(file=url)
if(!exists('DATA')){
    res = 'NODATA'
      stop('No answer from APPRIS Database : non-coding RNA requested or no
      access to the site')
        }
        json_data <- function(DATA)
        {
            res <- NULL
              for (i in 1:length(DATA[]))
                  {
                      res <- rbind(res, DATA[[i]][c('gene_id', 'end', 'start',
                      'transcript_id', 'type', 'annotation', 'reliability')])
                        }
                          return(res)
        }
        res <- as.data.frame(json_data(DATA))
        principal1 <- subset(res, stri_detect_fixed(res$reliability,
        'PRINCIPAL:1'))
        principal2 <- subset(res, stri_detect_fixed(res$reliability,
        'PRINCIPAL:2'))
        principal3 <- subset(res, stri_detect_fixed(res$reliability,
        'PRINCIPAL:3'))
        principal4 <- subset(res, stri_detect_fixed(res$reliability,
        'PRINCIPAL:4'))
        principal5 <- subset(res, stri_detect_fixed(res$reliability,
        'PRINCIPAL:5'))

if (nrow(principal1) != 0){
    res <- principal1
} else if (nrow(principal2) != 0) {
    res <- principal2
} else if (nrow(principal3) != 0) {
    res <- principal3
} else if (nrow(principal4) != 0) {
    res <- principal4
} else {
    res <- principal5
}

res$size <- as.integer(res$end) - as.integer(res$start) + 1
res <- subset(res, res$size == max(res$size))
res <- unique(res$transcript_id)
finish_time = Sys.time()
time = finish_time - start_time"
    
appris_res = @rget res time
res = String(res[1])
println("APPRIS result : $res")
println("APPRIS time : $time")
return(res)
end # end of APPRIS function

find_longer_variant = function(gene_name, dico_transcriptome)

    nb_variants = length(filter((k,v) -> startswith(k, "$gene_name:"), dico_transcriptome))
    println("nb variants : $nb_variants")
    variants_dico = filter((k,v) -> startswith(k, "$gene_name:"), dico_transcriptome)
    println("variants dico : $variants_dico")
    variants_lengths = Vector{Any}()
    for (k,v) in variants_dico
        println(k)
        println(v)
      push!(variants_lengths, length(v))
    end

    longer_variant = split(keys(filter((k,v) -> length(v) == maximum(variants_lengths),variants_dico)[1]))[2]
    
    println("$gene_name : longer variant = $longer_variant")
    return($longer_variant)
end






#split each sequence of fasta input file into individual fastas
# WARNING work only with ensembl fasta descriptions!
if unannotated_option == false
    FastaReader("$fastafile") do fr
    for (desc, seq) in fr
      desc_array = split(desc)
      #println( desc_array[1])
      gene_name = replace(desc_array[7], "gene_symbol:", "")
      ensembl_transcript_name = split(desc_array[1],'.')[1]
      ensembl_gene_name = split(replace(desc_array[4], "gene:",""), '.')[1]
      #println("ensembl gene name : $ensembl_gene_name")
      #println("ensembl transcript name : $ensembl_transcript_name")
      #println("gene name : $gene_name")
      #println("ensembl gene name : $(split(ensembl_gene_name, '.'))[1]")
      #println("test to see if gene name in list")
      if isempty(select_option) || (select_option != nothing && (gene_name in select_option || ensembl_transcript_name in select_option || ensembl_gene_name in select_option))         # take all sequence corresponding to asked gene names (select option)
        # take all if select_option not provided
        if(APPRIS_option != nothing) 
        APPRIS_transcript = APPRIS_function(ensembl_gene_name) 
        println("APPRIS selected variant : $APPRIS_transcript")
        end
      if APPRIS_option == nothing || (APPRIS_option != nothing && ensembl_transcript_name == APPRIS_transcript) || (APPRIS_option != nothing && APPRIS_transcript == "NODATA" && ensembl_transcript_name == find_longer_variant(gene_name,dico_transcriptome))
      println("$ensembl_transcript_name")

      if length("$seq") >= kmer_length 
        println("$gene_name:$ensembl_transcript_name : good length, continue")
        FastaWriter("$output/sequences/$kmer_length/$gene_name:$ensembl_transcript_name.fa") do fwsequence
        #for (desc2, seq2) in fr
        write(fwsequence, [">$gene_name:$ensembl_transcript_name", "$seq"])
        end
      
      else println("$gene_name:$ensembl_transcript_name : wrong length, ignored")
      end
    end
    end
    end
    end
  # if the input fasta is not a part of ensembl annotation (-u)
  else
    FastaReader("$fastafile") do fr
    for (desc, seq) in fr
    gene_name = desc
    if length("$seq") >= kmer_length
      FastaWriter("$output/sequences/$kmer_length/$desc.fa") do fwsequence
      write(fwsequence, [">$desc", "$seq"])
    end
  else println("$desc ========> wrong length, ignored")
    end
    end
end
end
println(readdir("$output/sequences/$kmer_length/"))
#println(keys(dico_transcriptome))

@passobj 1 workers() dico_transcriptome
@passobj 1 workers() unannotated_option
@passobj 1 workers() genome
@passobj 1 workers() transcriptome
@passobj 1 workers() level
@passobj 1 workers() output
@passobj 1 workers() kmer_length
@passobj 1 workers() stringent_option
@passobj 1 workers() admission_threshold


#Create the function f for extraction of specifics kmers for one sequence fasta
#file.

#print cleared lines used by replace_line function

for _ in 1:nworkers()
print("\n")
end

@everywhere function f(splitted_fasta_files)
#store worker id and X variable, for line replacement (if one core, one line up, if worker 3, 2 lines up)
if(myid() == 1) X = myid() else X = myid()-1 end
sleep(2)


enlapsed_time = 0
step_time = 0
kmers_analysed = 0


replace_line = function(nbline::Int64, x::String)
print(string("\u1b[$(nbline)F\u1b[2K"*"$x"*"\u1b[$(nbline)E"))
end

#println("$splitted_fasta_files")
#println(dico_transcriptome)

# manage the name of transcript/gene
  if unannotated_option == false # if annotated
    gene_name = split(splitted_fasta_files, ":")[1]
#    println("$gene_name")
    transcript_name = split(splitted_fasta_files, ":")[2]
#    test = filter((k,v) -> startswith(k, "$gene_name:"), dico_transcriptome)
#    println(test)
    nb_variants = length(filter((k,v) -> startswith(k, "$gene_name:"), dico_transcriptome))
#    println("nb variants : $nb_variants")
    variants_dico = filter((k,v) -> startswith(k, "$gene_name:"), dico_transcriptome)
    #println("variants dico : $variants_dico")
#    variants_lengths = Vector{Any}()
#    for (k,v) in variants_dico
#        println(k)
#        println(v)
#      push!(variants_lengths, length(v))
#    end

#    longer_variant = split(keys(filter((k,v) -> length(v) == maximum(variants_lengths),variants_dico)[1]))[2]
    
    #println("$gene_name : longer variant = $longer_variant")

  else # if unannotated
    gene_name = "$splitted_fasta_files"
    transcript_name = "$splitted_fasta_files"
  end
#println("level : $level")
if level == "gene"
  tag_file = "$gene_name-$transcript_name-gene_specific.fa"
end
if level == "transcript"
  tag_file = "$gene_name-$transcript_name-transcript_specific.fa"
end
if level == "chimera"
  tag_file = "$gene_name-chimera_specific.fa"
end


fasta_array = Array([])



# take the transcript sequence for jellyfish query
sequence_fasta = "error if seen"
FastaReader("$output/sequences/$kmer_length/$splitted_fasta_files") do fr
for (desc, seq) in fr
#    println("$desc")
  sequence_fasta = "$seq"
end
end
#println("$sequence_fasta")
#println("")

remotecall(replace_line, 1 , X, "worker number $(myid()), for $level $gene_name, reporting for duty ! :)\n")

#global genome
  kmercounts_genome = readstring(`jellyfish query -s "$output/sequences/$kmer_length/$splitted_fasta_files" "$genome"`)
remotecall(replace_line, 1 , X, "jellyfish query on genome.jf finished")
  kmercounts_transcriptome = readstring(`jellyfish query -s "$output/sequences/$kmer_length/$splitted_fasta_files" "$transcriptome"`)

remotecall(replace_line, 1 , X, "jellyfish query on transcriptome.jf finished")
#println("jellyfish query on transcriptome.jf finished")



kmercounts_genome = split(kmercounts_genome, "\n")[1:end-1]
kmercounts_genome_dico = Dict()
for mer in kmercounts_genome
mer = split(mer)
counts = mer[2]
seq = mer[1]
kmercounts_genome_dico["$seq"] = counts
end

kmercounts_transcriptome_dico = Dict()
kmercounts_transcriptome = split(kmercounts_transcriptome, "\n")[1:end-1]
for mer in kmercounts_transcriptome
mer = split(mer)
counts = mer[2]
seq = mer[1]
kmercounts_transcriptome_dico["$seq"] = counts
end



#initialisation of count variables
i = 0
total_kmers = length(keys(kmercounts_transcriptome_dico))

## create a new dictionary containing starts of kmers coordinates, to index
#kmers by their place on the sequence
kmer_starts = Dict()
kmer_placed = 0
#println("start processing kmer place on input sequence")
total_kmer_number = length(keys(kmercounts_transcriptome_dico))
for mer in keys(kmercounts_transcriptome_dico)
kmer_placed = kmer_placed +1

remotecall(replace_line, 1 , X, "kmers placed = $kmer_placed on $total_kmer_number")
kmer_interval = search("$sequence_fasta", "$mer")
kmer_start = minimum(collect(kmer_interval))
#kmer_end = maximum(collect(kmer_interval))
kmer_starts["$mer"] = kmer_start
end

#for mer in keys(kmercounts_transcriptome_dico)
for tuple in sort(collect(zip(values(kmer_starts),keys(kmer_starts))))
##for now, I can only have a sorted array of tuples, I have to extract sequence
#for each tuple
mer = last(tuple)
tic()
kmers_analysed = kmers_analysed + 1
per = round(kmers_analysed/total_kmers*100)
enlapsed_time = enlapsed_time + step_time
if enlapsed_time == 0
  time_remaining = 0
else
  time_remaining = Integer(round((total_kmers - kmers_analysed)/(kmers_analysed/enlapsed_time)))
end
ptime = Dates.canonicalize(Dates.CompoundPeriod(Dates.Second(time_remaining)))
#println("time : $ptime")  
remotecall(replace_line, 1 , X, "$level $gene_name : $kmers_analysed kmers analysed,$i specifics,( $per %) :remaining time -> $ptime")

#println("-------------------------------------------------------------------------------------------")

if !any(x->x == "$mer", keys(kmercounts_genome_dico))

    #println("$mer not defined : reversed ?")
  d = Dict("A"=>"T", "C"=> "G", "G"=>"C", "T"=>"A")
  new_mer = reverse(replace(mer, r"[ACGT]{1}", x->d[x]))  # "TGCA"
  #println("$(String(reverse_complement!(dna"$mer"))))")
  genome_count = kmercounts_genome_dico["$new_mer"]
  #println("genome count : $genome_count")
else

genome_count = kmercounts_genome_dico["$mer"]

end

transcriptome_count = kmercounts_transcriptome_dico["$mer"]
#println("transcriptome count : $transcriptome_count")
if level == "gene"
         #println("genome count : $genome_count")
          if genome_count == "1" || genome_count == "0"#if the kmer is present/unique or does not exist (splicing?) on the genome
            mer_regex = Regex(mer)
            variants_containing_this_kmer = keys(filter((k,v) -> ismatch(mer_regex, v), variants_dico))
            if stringent_option == true && float(transcriptome_count) == float(nb_variants) == float(length(variants_containing_this_kmer)) 
                #println("specific kmer found !")
              i = i+1
                tmp = length(variants_containing_this_kmer)
              push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)")
              push!(fasta_array,"$mer")

          elseif stringent_option == false && float(transcriptome_count) == float(length(variants_containing_this_kmer)) && float(transcriptome_count) > float(nb_variants*admission_threshold)

                #println("specific kmer found")
                i = i+1
                tmp = length(variants_containing_this_kmer)
              push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)" )
              push!(fasta_array, "$mer")

            elseif unannotated_option == true && float(transcriptome_count) == float(0)

              i = i+1
              push!(fasta_array,">$gene_name-$transcript_name.kmer$i" )
              push!(fasta_array,"$mer")

          end #end of fasta writing function
        
          #FastaWriter("$output/tags/$kmer_length/$tag_file") do fw
 #         writedlm("$output/tags/$kmer_length/$tag_file", transpose(fasta_array), "\n")
        end
  end


  #println("$genome_count")
if level == "transcript"
  #println("$genome_count")
  #println("$transcriptome_count")
  #println(parse(Int, genome_count) <= 1 )
  #println(float(transcriptome_count) == float(0))
  if unannotated_option == true && float(transcriptome_count) == float(0) && parse(Int, genome_count) <= 1
      i = i+1
              push!(fasta_array,">$gene_name.kmer$i")
              push!(fasta_array,"$mer")

  elseif float(transcriptome_count) == float(1) && parse(Int, genome_count) <= 1
    i = i+1
              push!(fasta_array,">$gene_name.kmer$i")
              push!(fasta_array,"$mer")
  end
 end

if level == "chimera" && unannotated_option == true
    #println("$genome_count")
    #println("$transcriptome_count")

    #println(parse(Int, genome_count) == 0 )
    #println(float(transcriptome_count) == float(0))
  if float(transcriptome_count) == float(0) && parse(Int, genome_count) == 0
    #println("YES")
    i = i+1
              push!(fasta_array,">$gene_name.kmer$i")
              push!(fasta_array,"$mer")
  end
  end


step_time = Integer(round(toq()))
#println("step time : $step_time")
end # end of kmer analysis
#println("$level $gene_name : $i specific kmers found")

writedlm("$output/tags/$kmer_length/$tag_file",
reshape(fasta_array,length(fasta_array)), "\n")
end #end of function
splitted_fasta_files = readdir("$output/sequences/$kmer_length/")

#pmap(f, splitted_fasta_files, unannotated_option, genome, transcriptome, level, output, kmer_length, stringent_option, admission_threshold)
pmap(f, splitted_fasta_files)

