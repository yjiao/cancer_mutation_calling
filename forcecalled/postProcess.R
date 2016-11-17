setwd("/seq/hacohenlab1/yjiao/prepost/voting/forcecalled")
files <- list.files(pattern=".*.maf")

si <- readRDS('/seq/hacohenlab1/yjiao/prepost/seqInfo/si.rds')

fh_info <- lapply(files, function(x) {
    print(x)
    tab <- read.table(x, header=T, sep='\t', quote='')
    if (nrow(tab) > 0) {
	x <- gsub(".maf", "", x)
	tab$Collaborator.Sample.ID <- x
	tab$build <- 37
	tumor <- si[si$Collaborator.Sample.ID==x, 'Broad.Sample.ID']
	tab$tumor_barcode <- tumor
	
	pid <- si[si$Collaborator.Sample.ID==x, 'ptid']
	idx <- (si$ptid == pid) & (si$order == 0)
	
	normal <- si[idx, 'Broad.Sample.ID']
	tab$normal_barcode <- normal
	
	write.table(tab, paste0("postprocessed/", x, "_postprocessed.maf"),
	    row.names=F, col.names=T, sep='\t', quote=F)

	pair <- paste0(pid, 'T-TP-NB-', tumor, '-', normal)
	fileloc <- paste0("/seq/hacohenlab1/yjiao/prepost/voting/forcecalled/postprocessed/", x, "_postprocessed.maf")
	tmp <- data.frame(pair_id = pair, maf_fc_vote=fileloc, individual_id=paste0(pid, 'T'))
    } else {
	print(paste("NO MUTATIONS FOUND:", x))
	tab$Collaborator.Sample.ID <- tab$chr
	tab$build <- tab$chr
	tab$tumor_barcode <- tab$chr
	tab$normal_barcode <- tab$chr
	write.table(tab, paste0("postprocessed/", x, "_postprocessed.maf"),
	    row.names=F, col.names=T, sep='\t', quote=F)
    }
})

fh_info <- do.call(rbind, fh_info)
write.table(fh_info, paste0("postprocessed/fh_info.tsv"),
	row.names=F, col.names=T, sep='\t', quote=F)

