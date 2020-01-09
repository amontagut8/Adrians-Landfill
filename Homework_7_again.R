####################Homework 7 redux

library(TBX20BamSubset)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(GenomicAlignments)
library(org.Mm.eg.db)
library(tools)
library(edgeR)

#Get all our data together and pretty it up

bam.files <- getBamFileList()
bam.files <- basename(bam.files)
bam.files <- file_path_sans_ext(bam.files)


#Assign each BAM file to an experimental condition (wild-type vs knockout)
experimental.conditions <- c("Wild-type", "Knockout")
experimental.conditions <- rep(experimental.conditions, each=3)

experimental.conditions <- factor(experimental.conditions)

names(experimental.conditions) <- bam.files

#Reference the genome 
mouse.genome <- TxDb.Mmusculus.UCSC.mm9.knownGene

#We can pull the levels of the seqname factor from this genome by using the seqlevels() function. 
#These refer to the levels of the seqname factor. So we need to reassign the "levels" vector to contain only our chrosome of interest

#We overwrrite this though, since we are only interested in chromosome 19
seqlevels(mouse.genome) <- c("chr19")

#Now for the fun stuff, get all the genes and transcripts from the genome file
chr19.genes <- genes(mouse.genome)
chr19.transcripts <- transcripts(mouse.genome)
chr19.exons <- exons(mouse.genome)

#Pull the gene ids

chr19.gene.ids <- chr19.genes$gene_id

#Start off our data frame as an nx1 df, with the Entrez IDs as our row names
count.table <- data.frame(row.names=chr19.gene.ids)


#Set the strandedness when putting together the proper GRanges object
strand.type <- "*" #Refers to both + and - strand. This is the default when setting up a GRanges object

bam.files.index <- c(1:length(bam.files))


#Generate the actual count table
bam.filepaths <- getBamFileList()


for (x in bam.files.index) {
  
  mouse.aligned.reads <- readGAlignments(bam.filepaths[x])
  mouse.aligned.reads.ranges <- ranges(mouse.aligned.reads)
  #Generate a GRanges object from each BAM file and then count the overlaps with the annotation file
  
  mouse.aligned.reads.object <- GRanges(seqnames = "chr19", ranges = mouse.aligned.reads.ranges, strand = strand.type)
  mouse.overlaps <- countOverlaps(chr19.genes, mouse.aligned.reads.object)
  mouse.overlaps.data.frame <- data.frame(mouse.overlaps)
  
  names(mouse.overlaps.data.frame) <- experimental.conditions[x]
  
  count.table <- cbind(count.table, mouse.overlaps.data.frame)
  
}

#Prune the actual count table to remove genes where no transcript was detected at all across all runs
#We can find this by removing rows where the sum of all overlaps is 0
count.table <- cbind(count.table, 'Transcript_Sums'=rowSums(count.table))

count.table.pruned <- count.table[count.table$Transcript_Sums != 0,]

#So now we need an annotation file for converting Entrez ID names into HGNC gene symbols
annotation <- select(org.Mm.eg.db, rownames(count.table.pruned), keytype = "ENTREZID", "SYMBOL")

annotation.geneIDs <- annotation$ENTREZID
annotation.geneSymbols <- annotation$SYMBOL

#Are there any NaNs? 

any(is.na(mapping.annotation.geneSymbols))

#Returns TRUE. So we need to remove these rows.
annotation <- annotation[is.na(annotation$SYMBOL)==FALSE,]

#And now we finish building the count table by adding the gene symbol column to the count table
count.table.pruned["ENTREZ_ID"]=row.names(count.table.pruned)


count.table.final <- merge(x=count.table.pruned, y = annotation,
                           by.x="ENTREZ_ID", by.y="ENTREZID")



#So now we have the count table, we can figure out which genes were most differentially regulated between wild-type and knockouts


#But first we need to get it into a format that we can feed it into a DGEList object

count.table.edge <- count.table.final[2:7]
row.names(count.table.edge) <- count.table.final$SYMBOL


#Fun stuff to do with a DGEList
gene.expression <- DGEList(counts = count.table.edge, group = experimental.conditions)
gene.expression <- calcNormFactors(gene.expression)

question.2.gene.expression <- estimateCommonDisp(gene.expression)
question.2.gene.expression <- estimateTagwiseDisp(gene.expression)


#The meat of the matter - how to do a DGA
differential.expression <- exactTest(question.2.gene.expression)

print(differential.expression)

#Now to get the top five genes, we can go about this in different ways. I would have gone for the 5 lowest p-values
#But there is another function topTags

differential.expression.top <- topTags(differential.expression, n=5, sort.by="PValue")

print("The top five most differentially expressed genes are:")
print(row.names(differential.expression.top))










