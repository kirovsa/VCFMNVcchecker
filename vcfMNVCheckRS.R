#Evaluate vcf file for MNVs
library(shiny)
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)
#install.packages("shinyWidgets")
library(shinyWidgets)
library(shinyjs)
options(shiny.maxRequestSize = 100*1024^2)


checkVCF<-function(x,vcf) { #x is a missense row from 
  check<-1
  checkReason<-""
  #are there any rows with missense?
  if (length(x)==0) {
    check<-0
    checkReason<-"No parseable annotation"
    print(checkReason)
    return(c(check,checkReason))
  }
  #Check if annotation is present, specifically substtitution in protein coordinates
  if (length(grep("p\\.[\\w\\d]+",x[8],perl=T))==0) {
    #No annotation
    check<-0
    checkReason<-"No mutation annotation"
    print(checkReason)
    return(c(check,checkReason))
  }
  #Check if allele frequency can be parsed
  if (sum(vcf[['af']],vcf[['fa']])==0 & length(vcf[['strelka']])==0) {
    check<-0
    checkReason<-"No parseable allele frequency data"
    print(checkReason)
    return(c(check,checkReason))
  }
  return(c(check,checkReason))
}

bindRows<-function (resmatrix,inputrow) {
  resmatrix<-rbind(resmatrix,inputrow)
}

getSubstitution  <- function (infostr) {
  realsub<-"---"
  
  if (length(grep("p\\.\\w{3,3}\\d+\\w{3,3}",infostr,perl=T))>0) {
    realsub<-str_extract(infostr[grep("p\\.\\w{3,3}\\d+\\w{3,3}",infostr,perl=T)],"p\\.\\w{3,3}\\d+\\w{3,3}")
  }
  if (length(realsub)==0) {
    realsub<-"*"
  }
  return(realsub)
}

#checkRes<-checkVCF(vcf[['vcfpass']][1,],)

getFrequency<-function(rowD,vcf) {
  #Get the allele frequency per sample
  freq.df<-data.frame(sample=c(),frequency=c(),count=c())
  snames<-vcf[['samples']]
  for (i in 10:(9+length(snames))) {
    alleleD<-strsplit(as.character(rowD[i]),"\\:")
    if (length(vcf[['dp']])!=0) { 
      dpD<-as.numeric(alleleD[[1]][vcf[['dp']]])
    }
    if (length(vcf[['ad']])!=0) {
      adD<-strsplit(alleleD[[1]][vcf[['ad']]],"\\,")[[1]]
      dpD<-as.numeric(adD[1])+as.numeric(adD[2])
    }
    if (length(vcf[['strelka']])>1) {
      refBase<-paste0(rowD[4],"U")
      altBase<-paste0(rowD[5],"U")
      altind<-grep(altBase,vcf[['format']])
      refind<-grep(refBase,vcf[['format']])
      #Tier1 only
      refDat<-as.numeric(strsplit(alleleD[[1]][refind],"\\,")[[1]][1])
      altDat<-as.numeric(strsplit(alleleD[[1]][altind],"\\,")[[1]][1])
      saf<-altDat/(altDat+refDat)
      df<-data.frame(sample=i,frequency=saf,count=(altDat+refDat))
      freq.df<-rbind(freq.df,df)
    }
    if (length(vcf[['af']])!=0 | length(vcf[['fa']])!=0) {
      saf<-as.numeric(alleleD[[1]][sum(vcf[['af']],vcf[['fa']], na.rm = TRUE)])
      df<-data.frame(sample=i,frequency=as.numeric(saf),count=dpD)
      freq.df<-rbind(freq.df,df)
    }
    if (length(vcf[['ad']])!=0) {
      
    }
    if (length(vcf[['ac']])!=0) {
      
    }
  }
  return(freq.df)
}

getTranscriptCodon<-function (x,vcf) {
  
  #Sample data
  AFperS<-getFrequency(x,vcf)
  
  #Transcript data
  transcriptD<-strsplit(tolower(x[8]),"\\|transcript\\|")
  if (length(transcriptD)==0) return(NA)
  cold<-lapply(transcriptD,function(y) strsplit(y,"\\|"))
  gd<-strsplit(cold[[1]][[1]],"\\|")
  gene<-gd[[length(gd)]]
  
  cold[[1]][[1]]<-NULL
  perT<-lapply(cold[[1]],function(y) y[1])
  perC<-lapply(cold[[1]],getSubstitution )
  if (length(perT)==0) {
    print(gene)
    return(NULL)
    
  }
  
  if (length(as.vector(unlist(perT)))!=length(as.vector(unlist(perC)))) {
    print(gene)
    return(NULL)
  }
  res<-data.frame(Transcript=as.vector(unlist(perT)),substitution=as.vector(unlist(perC)))
  
  position<-str_extract(res$substitution,"\\d+")
  res$position<-position
  res$Gene<-gene
  res<-merge(res,AFperS,all=TRUE)
  return(res)
}


combineMNVs<-function(rowsMNV) {
  #Assume 2 only for now
  #If these are the same substitutions they may be merged or otherwise ineligible
  if (rowsMNV[1,8]==rowsMNV[2,8]) {
    return(NULL)
  }
  rowData<-c(as.character(rowsMNV[1,c(1,2,3,4,5,8,12,14)]),as.character(rowsMNV[2,c(8,12,14)]))
  names(rowData)<-c("UID","Gene","Sample","Transcript","Position","Mutation1","Freq1","Reads1","Mutation2","Freq2","Reads2")
  return(rowData)
}

getMNVs<-function(uidlist,allcodons.mnv) {
  mnv.data<-lapply(uidlist, function(x) combineMNVs(allcodons.mnv[allcodons.mnv$UID==x,]))
  return(Filter(Negate(is.null), mnv.data))
}

readVCF<-function(inputvcf) {
  tmp_vcf<-readLines(inputvcf,skipNul = T)
  mnvdetected<-grep("\\;type=mnv\\;",tmp_vcf,ignore.case = T)
  if (length(mnvdetected)>0) {
    #We have MNVS, we are done!
    vcf<-list(check=c(3,"MNVs are correctly detected"),strelka=NA,samples=NA,format=NA,af=NA,ad=NA,ac=NA,fa=NA,dp=NA,isfiltered=NA,vcf=NA,vcfpass=NA)
    return(vcf)
  }
  tmp_vcf_data<-read.table(inputvcf, stringsAsFactors = FALSE,quote = "",skipNul = T)
  #Do we have strelka?
  strelka<-grep("strelka",tmp_vcf)
  # filter for the columns names
  tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
  vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
  names(tmp_vcf_data)<-vcf_names
  vcf<-tmp_vcf_data[grep("missense|synonymous",tmp_vcf_data[,8]),]
  vcf.pass<-vcf[vcf[,7]=="PASS",]
  snames<-vcf_names[10:length(vcf_names)]
  formatD<-unlist(strsplit(vcf.pass[1,9],"\\:"))
  af<-grep("\\bAF\\b",formatD)
  ad<-grep("\\bAD\\b",formatD)
  ac<-grep("\\bAC\\b",formatD)
  fa<-grep("\\bFA\\b",formatD)
  dp<-grep("\\bDP\\b",formatD)
  isfiltered<-dim(vcf.pass)[1]==dim(vcf)[1]
  vcf<-list(strelka=strelka,samples=snames,format=formatD,af=af,ad=ad,ac=ac,fa=fa,dp=dp,isfiltered=isfiltered,vcf=vcf,vcfpass=vcf.pass)
  checkRes<-checkVCF(vcf[['vcf']][1,],vcf)
  vcf[['check']]<-checkRes
  return(vcf)
}


analyzeVCF<-function (VCFfile,pass) {
  vcf<-readVCF(VCFfile)
  #We failed to parse, missing info, etc.
  if (vcf[['check']][[1]]!=1) {
    return(vcf)
  }
  print(paste("File read",dim(vcf[['vcf']])))
  #Will be option- vcf pass or vcf
  if (pass) {
    vcfD<-vcf[['vcfpass']]
  }
  else {
    vcfD<-vcf[['vcf']]
  }
  allcodons<-apply(vcfD,1,function(x) getTranscriptCodon(x,vcf))
  allcodons.df<-bind_rows(allcodons)
  allcodons.df<-allcodons.df[allcodons.df$substitution!="---",]
  allcodons.df$Reads<-allcodons.df$frequency*allcodons.df$count
  allcodons.df<-allcodons.df[allcodons.df$Reads>3,]
  allcodons.df$UID=interaction(allcodons.df$Transcript,allcodons.df$sample,allcodons.df$position)
  allcodons.df$UID<-as.character(allcodons.df$UID)
  print(paste("Data frame ready",dim(allcodons.df)))
  allcodons.stats<-aggregate(substitution ~ Gene + sample + Transcript + position,data=allcodons.df,length)
  allcodons.stats$UID=interaction(allcodons.stats$Transcript,allcodons.stats$sample,allcodons.stats$position)
  allcodons.stats$UID<-as.character(allcodons.stats$UID)
  allcodons.mnv<-merge(allcodons.stats[allcodons.stats$substitution==2,],allcodons.df,by="UID")
  #mnvs.freq<-dcast(allcodons.mnv, Transcript.x + sample.x + position.x  ~ substitution.y, value.var="frequency")
  mnvs<-getMNVs(unique(as.character(allcodons.mnv$UID)),allcodons.mnv)
  #mnvs.df<-bind_rows(mnvs)
  mnvs.df<-as.data.frame(do.call(rbind, mnvs))
  colnames(mnvs.df)<-c("UID","Gene","Sample","Transcript","Position","Mutation1","Freq1","Reads1","Mutation2","Freq2","Reads2")
  #We have to take 1 per gene! top1 should do it.
  mnvs.df$GUID<-as.character(interaction(mnvs.df$Gene,mnvs.df$Sample,mnvs.df$Mutation2,mnvs.df$Position))
  mnvs.df<-mnvs.df[match(unique(mnvs.df$GUID), mnvs.df$GUID),]
  mnvs.df$Reads<-unlist(apply(mnvs.df[,c("Reads1","Reads2")], 1, FUN=min))
  vcf[['mnv']]=mnvs.df
  vcf[['dims']]<-c(dim(allcodons.df),dim(allcodons.stats),dim(allcodons.mnv),dim(mnvs.df))
  return(vcf)
}

ui <- fluidPage(
  useShinyjs(),
  # App title ----
  titlePanel("VCF MNV test applicaction"),
  p("Explore if incorrectly annotated MNVs are present in a VCF"), 
  p("Please wait until graph is rendered or error message appears. If none of this happens in more than 10 minutes please send us an example of your file"),
  # Sidebar layout with input and output definitions ----
  fluidRow(
    column(10,
           tableOutput('firstMNVs')
    )
  ),
  
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput(
        "VCFfile",
        "VCF",
        multiple = FALSE,
        accept = ".vcf",
        buttonLabel = "Browse...",
        placeholder = "No file selected"
      ),
      
      # Sample selector
      selectInput(
        inputId = "Sample",
        label = "Sample index",
        NULL,
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      switchInput(
        "passSwitch",
        label = "Apply pass filter",
        value = FALSE,
        onLabel = "PASS",
        offLabel = "ALL"
      )
      
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: VCF MNV plot ----
      plotOutput(outputId = "VCFplot"),
      textOutput(outputId = "messages"),
      textOutput(outputId = "VCFcharacteristics")
    )
  )
)

vcfFileUpload<-function(session,input){
  file <- input$VCFfile
  print(file$datapath)
  ext <- tools::file_ext(file$datapath)
  passFlag<-input$passSwitch
  req(file) 
  validate(need(ext == "vcf", "Please upload a VCF file"))
  mnv<-analyzeVCF(file$datapath,passFlag)
  #if reading fails
  if (mnv[['check']][1]!=1){
    return(mnv)
  }
  samples<-c("ALL",mnv[["samples"]])
  print(paste("SAMPLES***",samples))
  updateSelectInput(session, "Sample",  choices = samples,
                    selected = samples[1])
  return(mnv)
}

plotMNVObject <-function(mnv,input) {
  print(paste("In render",mnv[['dims']]))
  #Samples is index vs name?
  print(paste("Samples",unique(mnv[['mnv']]$Sample)))
  sample<-grep(input$Sample,mnv[["samples"]])+9
  print(paste("For sample",input$Sample,"with index",sample))
  if (length(sample)==0) {
    mnv.dat.plot<-mnv[['mnv']]
  }
  else {#All samples
      mnv.dat.plot<-mnv[['mnv']][mnv[['mnv']]$Sample==sample,]
  }
  if(is.null(mnv.dat.plot)) {
    print("MNV object is null, returning")
    return(list(rplot=NA,rcor=NA,scor=NA,error="No file uploaded yet"))
  }
  print("We have the object")
  print(paste("Potential MNVS dataframe",dim(mnv.dat.plot),"for sample",input$Sample))
  if(dim(mnv.dat.plot)[1]<3) {
    print("Not enough SNVs in the same codon")
    return(list(rplot=NA,rcor=NA,scor=NA,error="Not enough SNPs in the same codon to determine MNV impact",mnv=mnv.dat.plot))
  }
  else {
    print(mnv.dat.plot)
    mnv.dat.plot$Freq1<-as.numeric(as.character(mnv.dat.plot$Freq1))
    mnv.dat.plot$Freq2<-as.numeric(as.character(mnv.dat.plot$Freq2))
    mnv.dat.plot$Reads<-as.numeric(as.character(mnv.dat.plot$Reads))
    pcor<-round(cor(mnv.dat.plot$Freq1,mnv.dat.plot$Freq2),2)
    scor<-round(cor(mnv.dat.plot$Freq1,mnv.dat.plot$Freq2,method="spearman"),2)
    if (is.na(pcor)) {
      fres<-list(rplot=NA,pcor=pcor,scor=scor,error="Possible redundant data, correlation cannot be computed")
      return(fres)
    } 
    print(paste("Pearson",pcor))
    print(paste("Spearman",scor))
    
    rplot<-ggplot(mnv.dat.plot,aes(x=Freq1,y=Freq2,col=Reads))+geom_point() + ggtitle(paste("Same codon SNVs for sample",input$Sample)) +
      annotate(geom="text", x=0.2, y=0.6, label=paste("Spearman",scor),color="red") +
      annotate(geom="text", x=0.2, y=0.63, label=paste("Pearson",pcor),color="red") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
    fres<-list(rplot=rplot,pcor=pcor,scor=scor,error="",mnv=mnv.dat.plot)
    return(fres)
  }
}


#Plot per sample and per PASS/no filter if possible

server <- function (input, output, session) {
  
  rmnv<-reactiveValues(mnv=NULL,tracker=1)
  #only on file upload
  
  observeEvent(input$VCFfile, {
    print("Calculating")
    hide("VCFplot")
    tmnv<-vcfFileUpload(session,input)
    print(tmnv[['check']][1])
    if (tmnv[['check']][1]!=1) {
      print(tmnv[['check']][2])
      output$messages<-renderText(tmnv[['check']][2])
      return(tmnv)
    }
    print(paste("In observer:",tmnv[['mnv']]))
    rmnv$mnv <- tmnv 
    rmnv$tracker=2
  })
  
  toListen <- reactive({
    list(input$Sample,input$passSwitch)
  })
  
  observeEvent(toListen(), {  
    fres<-plotMNVObject(rmnv$mnv,input)
    if (fres[['error']]=='') { 
      if (!is.na(fres[['rplot']])) {
        print("Printing the plot")
        output$VCFplot<-renderPlot(fres[['rplot']])
        output$firstMNVs<-renderTable(fres[['mnv']])
        show("VCFplot")
      }
      print("Checking correlations")
      if (fres[['pcor']]>0.9 & fres[['scor']]>0.9) {
        decision<-"VCF contains incorrectly annotated MNVs"
      }
      else {
        decision<-"VCF does not contain incorrectly annotated MNVs"
      }
        print("Done with decision")
        textRes<-paste(decision,"Pearson correlation is ",fres[['pcor']],"and Spearman is",fres[['scor']])
        output$messages<-renderText(textRes)
      }
    else {
      output$messages<-renderText(fres[['error']])
      output$firstMNVs<-renderTable(fres[['mnv']])
      hide("VCFplot")
    }
  })
  
}


shinyApp(ui, server)

