db <- read.table("Results.txt",header=T,sep="\t")
names(db)

#ACHE
##Select the ACHE data
db_ache <- db[db$Gene=="ache",]
##Find rs variant from the list
intersect(db_ache$name,db$control_snp)
##Distance rs variants to start of the ACHE gene (100889994 )
ord_ache <- (db_ache$chromStart)-100889994
db_ache <- data.frame(db_ache,ord_ache)


#ALDOA
##Select the ALDOA data
db_aldoa <- db[db$Gene=="aldoa",]
##Find rs variant from the list
intersect(db_aldoa$name,db$control_snp)
##Distance rs variants closer to start of the ALDOA gene (30064274)
ord_aldoa <- (db_aldoa$chromStart)-30064274
db_aldoa <- data.frame(db_aldoa,ord_aldoa)


#ANAPC4
##Select the ANAPC4 data
db_anapc4 <- db[db$Gene=="anapc4",]
##Find rs variant from the list
intersect(db_anapc4$name,db$control_snp)
##Distance rs variants closer to start of the ANAPC4 gene (25377263)
ord_anapc4 <- (db_anapc4$chromStart)-25377263
db_anapc4 <- data.frame(db_anapc4,ord_anapc4)


#BCL7A
##Select the BCL7A data
db_bcl7a <- db[db$Gene=="bcl7a",]
##Find rs variant from the list
intersect(db_bcl7a$name,db$control_snp)
##Distance rs variants closer to start of the bcl7a gene (122021884)
ord_bcl7a <- (db_bcl7a$chromStart)-122021884

#BMP2
##Select the BMP2 data
db_bmp2 <- db[db$Gene=="bmp2",]
##Find rs variant from the list
intersect(db_bmp2$name,db$control_snp)
##Distance rs variants closer to start of the bmp2 gene (6767686)
ord_bmp2 <- (db_bmp2$chromStart)-6767686


#BPTF
##Select the BPTF data
db_bptf <- db[db$Gene=="bptf",]
##Find rs variant from the list
intersect(db_bptf$name,db$control_snp)
##Distance rs variants closer to start of the bpdf gene (67826033)
ord_bptf <- (db_bptf$chromStart)-67826033


#C6ORF106
##Select the C6ORF106 data
db_c6orf106 <- db[db$Gene=="c6orf106",]
##Find rs variant from the list
intersect(db_c6orf106$name,db$control_snp)
##Distance rs variants closer to start of the c6orf106 gene (34587288)
ord_c6orf106 <- (db_c6orf106$chromStart)-34587288
getwd()
write.table(ord_c6orf106,"clipboard-16384",sep="\t",col.names=NA)

#CASC20
##Select the CASC20 data
db_casc20 <- db[db$Gene=="casc20",]
##Find rs variant from the list
intersect(db_casc20$name,db$control_snp)
##Distance rs variants closer to start of the CASC20 gene (6446723)
ord_casc20 <- (db_casc20$chromStart)-6446723


#CYB5B
##Select the CASC20 data
db_cyb5b <- db[db$Gene=="cyb5b",]
##Find rs variant from the list
intersect(db_cyb5b$name,db$control_snp)
##Distance rs variants closer to start of the CYB5B gene (69424619)
ord_cyb5b <- (db_cyb5b$chromStart)-69424619


#DNAJC27
##Select the DNAJC27 data
db_dnajc27 <- db[db$Gene=="dnajc27",]
##Find rs variant from the list
intersect(db_dnajc27$name,db$control_snp)
##Distance rs variants closer to start of the DNAJC27 gene (24943642)
ord_dnajc27 <- (db_dnajc27$chromStart)-24943642


# ECE2
##Select the ECE2 data
db_ece2 <- db[db$Gene=="ece2",]
##Find intersection with list
intersect(db_ece2$name,db$control_snp)
##Distance rs variants closer to start of the ECE2 gene (184276022)
ord_ece2 <- (db_ece2$chromStart)-184276022
getwd()
write.table(ord_ece2,"clipboard-16384",sep="\t",col.names=NA)

#GIPC2
##Select the GIPC2 data
db_gipc2 <- db[db$Gene=="gipc2",]
str(db_gipc2)
##Find intersection with list
intersect(db_gipc2$name,db$control_snp)
##Distance rs variants closer to start of the GIPC2 gene (78045966)
ord_gipc2 <- (db_gipc2$chromStart)-78045966

#IGF2BP2
##Select the IGF2BP2 data
db_igf2bp2 <- db[db$Gene=="igf2bp2",]
##Find intersection with list
intersect(db_igf2bp2$name,db$control_snp)
##Distance rs variants closer to start of the IGF2BP2 gene (185643130)
ord_igf2bp2 <- (db_igf2bp2$chromStart)-185643130


#LOXL4
##Select the LOXL4 data
db_loxl4 <- db[db$Gene=="loxl4",]
##Find intersection with list
intersect(db_loxl4$name,db$control_snp)
##Distance rs variants closer to start of the LOXL4 gene (98247690)
ord_loxl4 <- (db_loxl4$chromStart)-98247690


#MC4R
##Select MC4R data
db_mc4r <- db[db$Gene=="mc4r",]
##Find intersection with list
intersect(db_mc4r$name,db$control_snp)
##Distance rs variants closer to start of the MC4R gene (60371062)
ord_mc4r <- (db_mc4r$chromStart)-60371062


#MIR1538
##Select MIR1538 data 
db_mir1538 <- db[db$Gene=="mir1538",]
##Find intersection with list
intersect(db_mir1538$name,db$control_snp)
##Distance rs variants closer to start of the MIR1538 gene (69565808)
ord_mir1538 <- (db_mir1538$chromStart)-69565808

#MTCH2
##Select MTCH2 data 
db_mtch2 <- db[db$Gene=="mtch2",]
##Find intersection with list
intersect(db_mtch2$name,db$control_snp)
##Distance rs variants closer to start of the MTCH2 gene (47617317)
ord_mtch2 <- (db_mtch2$chromStart)-47617317


#NCAN
##Select NCAN data 
db_ncan <- db[db$Gene=="ncan",]
##Find intersection with list
intersect(db_ncan$name,db$control_snp)
##Distance rs variants closer to start of the NCAN gene (19211958)
ord_ncan <- (db_ncan$chromStart)-19211958


NRXN3
#Select NRXN3 data
db_nrxn3 <- db[db$Gene=="nrxn3",]
##Find intersection with list
intersect(db_nrxn3$name,db$control_snp)
which(db_nrxn3$name=="rs10146997")
##Distance rs variants closer to start of the NRXN3 gene (78170373)
ord_nrxn3 <- (db_nrxn3$chromStart)-78170373


#PCCB
##Select PCCB data
db_pccb <- db[db$Gene=="pccb",]
##Find intersection with list
intersect(db_pccb$name,db$control_snp)
##Distance rs variants closer to start of the PCCB gene (136250355)
ord_pccb <- (db_pccb$chromStart)-136250355
getwd()
write.table(ord_pccb,"clipboard-16384",sep="\t",col.names=NA)

#PCSK1
##Select PCSK1 data
db_pcsk1 <- db[db$Gene=="pcsk1",]
##Find intersection with list
intersect(db_pcsk1$name,db$control_snp)
##Distance rs variants closer to start of the PCSK1 gene (96390333)
ord_pcsk1 <- (db_pcsk1$chromStart)-96390333

#SLC39A8
##Select PCSK1 data
db_slc39a8 <- db[db$Gene=="slc39a8",]
##Find intersection with list
intersect(db_slc39a8$name,db$control_snp)
##Distance rs variants closer to start of the SLC39A8 gene (102261677)
ord_slc39a8 <- (db_slc39a8$chromStart)-102261677
getwd()
write.table(ord_slc39a8,"clipboard-16384",sep="\t",col.names=NA)

#ZBTB7B
##Select ZBTB7B data
db_zbtb7b <- db[db$Gene=="zbtb7b",]
##Find intersection with list
intersect(db_zbtb7b$name,db$control_snp)
##Distance rs variants closer to start of the ZBTB7B gene (155002657)
ord_zbtb7b <- (db_zbtb7b$chromStart)-155002657
getwd()
write.table(ord_zbtb7b,"clipboard-16384",sep="\t",col.names=NA)





























































db_bcl7a <- data.frame(db_bcl7a,ord_bcl7a)
getwd()
write.table(db_anapc4,"clipboard-16384",sep="\t",col.names=NA)



















