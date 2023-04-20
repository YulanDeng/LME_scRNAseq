## Script to Build a neural network–based classifier for TCR-antigen specificity
##Yulan Deng, last updated 2023-4-20
##my e-mail:kndeajs@163.com

###################################
#basic process of training dataset#
###################################
#${smp} is the sample ID
#${workDir} is the working directory
#${analysis_name} is the name for analysis
#${clonotypeDir} is the directory for clonotypes from MixCR

#convert sra file to fastq file
fastq-dump --split-3 ${smp}

#MixCR pipeline of TCRseq
docker run -it --rm -m 4g -v ${workDir}:/work milaboratory/mixcr:latest \
mixcr analyze amplicon --species HomoSapiens --starting-material dna --5-end v-primers --3-end j-primers --adapters no-adapters --receptor-type TRB -t 8 ${smp}_1.fastq ${smp}_2.fastq ${analysis_name}

#positive and negative datasets
#Set working directory
setwd(workDir)


label <- c("pre","mo1","mo2","mo3")
tcr_res_list <- lapply(label,function(x){
	res <- read.table(file = paste0(clonotypeDir,x,"/",x,".clonotypes.TRB.txt"),
	sep="\t",quote = "",stringsAsFactors=F,header=T)
	reslabel <- paste(res[,"allVHitsWithScore"],res[,"allDHitsWithScore"],
	res[,"allJHitsWithScore"],res[,"nSeqCDR3"],res[,"aaSeqCDR3"],sep="@")
	resclonecount <- res[,"cloneCount"]
	resclonefraction <- res[,"cloneFraction"]
	res_df <- data.frame(reslabel=reslabel,resclonecount=resclonecount,
	resclonefraction=resclonefraction,stringsAsFactors=F)
	return(res_df)
})

#pretreatment and 2 month
inter13 <- intersect(tcr_res_list[[1]][,1],tcr_res_list[[3]][,1])
fre_threshold_down <- quantile(tcr_res_list[[1]][,3],probs=0.6)
neg <- tcr_res_list[[1]][(!(tcr_res_list[[1]][,1]%in%inter13))&(tcr_res_list[[1]][,2]>=2)&(tcr_res_list[[1]][,3] < fre_threshold_down),]
neg <- neg[1:3217,]
fre_threshold <- quantile(tcr_res_list[[3]][,3],probs=0.95)
pos <- tcr_res_list[[3]][(!(tcr_res_list[[3]][,1]%in%inter13))&(tcr_res_list[[3]][,2]>=10)&(tcr_res_list[[3]][,3]>=fre_threshold),]
pos_trb_cdr3 <- sapply(strsplit(pos[,1],"@"),function(x) x[5])
pos_trb_v_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[1],"*",fixed=T)[[1]][1])
pos_trb_d_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[2],"*",fixed=T)[[1]][1])
pos_trb_j_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[3],"*",fixed=T)[[1]][1])
pos_df <- data.frame(TRB_cdr3=pos_trb_cdr3,TRB_v_gene=pos_trb_v_gene,
TRB_d_gene=pos_trb_d_gene,TRB_j_gene=pos_trb_j_gene,Clonality="na",
pmhc_code="EGFRpL858R",id="KITDFGRAK",binds=1,TCR_id=seq(nrow(pos)),
stringsAsFactors=F)
pos_df <- pos_df[grep("_", pos_df[,1], invert = T,fixed=T),]
pos_df <- pos_df[grep("*", pos_df[,1], invert = T,fixed=T),]
neg_trb_cdr3 <- sapply(strsplit(neg[,1],"@"),function(x) x[5])
neg_trb_v_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[1],"*",fixed=T)[[1]][1])
neg_trb_d_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[2],"*",fixed=T)[[1]][1])
neg_trb_j_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[3],"*",fixed=T)[[1]][1])
neg_df <- data.frame(TRB_cdr3=neg_trb_cdr3,TRB_v_gene=neg_trb_v_gene,
TRB_d_gene=neg_trb_d_gene,TRB_j_gene=neg_trb_j_gene,Clonality="na",
pmhc_code="EGFRpL858R",id="KITDFGRAK",binds=0,TCR_id=seq(nrow(neg))+nrow(pos),
stringsAsFactors=F)
neg_df <- neg_df[grep("_", neg_df[,1], invert = T,fixed=T),]
neg_df <- neg_df[grep("*", neg_df[,1], invert = T,fixed=T),]
tran_df <- rbind(pos_df,neg_df)
write.table(tran_df, file="EGFRpL858R.tran_df.csv",sep=",", quote=FALSE, row.names=F, col.names=T)

#1 month
inter12 <- intersect(tcr_res_list[[1]][,1],tcr_res_list[[2]][,1])
fre_threshold_down <- quantile(tcr_res_list[[1]][,3],probs=0.6)
neg <- tcr_res_list[[1]][(!(tcr_res_list[[1]][,1]%in%inter12))&(tcr_res_list[[1]][,2]>=2)&(tcr_res_list[[1]][,3] < fre_threshold_down),]
neg <- neg[3218:4825,]
fre_threshold <- quantile(tcr_res_list[[2]][,3],probs=0.95)
pos <- tcr_res_list[[2]][(!(tcr_res_list[[2]][,1]%in%c(inter12,tcr_res_list[[3]][,1])))&(tcr_res_list[[2]][,2]>=10)&(tcr_res_list[[2]][,3]>=fre_threshold),]
pos_trb_cdr3 <- sapply(strsplit(pos[,1],"@"),function(x) x[5])
pos_trb_v_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[1],"*",fixed=T)[[1]][1])
pos_trb_d_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[2],"*",fixed=T)[[1]][1])
pos_trb_j_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[3],"*",fixed=T)[[1]][1])
pos_df <- data.frame(TRB_cdr3=pos_trb_cdr3,TRB_v_gene=pos_trb_v_gene,
TRB_d_gene=pos_trb_d_gene,TRB_j_gene=pos_trb_j_gene,Clonality="na",
pmhc_code="EGFRpL858R",id="KITDFGRAK",binds=1,TCR_id=seq(nrow(pos)),
stringsAsFactors=F)
pos_df <- pos_df[grep("_", pos_df[,1], invert = T,fixed=T),]
pos_df <- pos_df[grep("*", pos_df[,1], invert = T,fixed=T),]
neg_trb_cdr3 <- sapply(strsplit(neg[,1],"@"),function(x) x[5])
neg_trb_v_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[1],"*",fixed=T)[[1]][1])
neg_trb_d_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[2],"*",fixed=T)[[1]][1])
neg_trb_j_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[3],"*",fixed=T)[[1]][1])
neg_df <- data.frame(TRB_cdr3=neg_trb_cdr3,TRB_v_gene=neg_trb_v_gene,
TRB_d_gene=neg_trb_d_gene,TRB_j_gene=neg_trb_j_gene,Clonality="na",
pmhc_code="EGFRpL858R",id="KITDFGRAK",binds=0,TCR_id=seq(nrow(neg))+nrow(pos),
stringsAsFactors=F)
neg_df <- neg_df[grep("_", neg_df[,1], invert = T,fixed=T),]
neg_df <- neg_df[grep("*", neg_df[,1], invert = T,fixed=T),]
tran_df <- rbind(pos_df,neg_df)
write.table(tran_df, file="EGFRpL858R.tran_df.mo1.csv",sep=",", quote=FALSE, row.names=F, col.names=T)

#3 month
inter14 <- intersect(tcr_res_list[[1]][,1],tcr_res_list[[4]][,1])
fre_threshold_down <- quantile(tcr_res_list[[1]][,3],probs=0.6)
neg <- tcr_res_list[[1]][(!(tcr_res_list[[1]][,1]%in%inter14))&(tcr_res_list[[1]][,2]>=2)&(tcr_res_list[[1]][,3] < fre_threshold_down),]
neg <- neg[4826:6435,]
fre_threshold <- quantile(tcr_res_list[[4]][,3],probs=0.95)
pos <- tcr_res_list[[4]][(!(tcr_res_list[[4]][,1]%in%c(inter14,tcr_res_list[[3]][,1])))&(tcr_res_list[[4]][,2]>=10)&(tcr_res_list[[4]][,3]>=fre_threshold),]
pos_trb_cdr3 <- sapply(strsplit(pos[,1],"@"),function(x) x[5])
pos_trb_v_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[1],"*",fixed=T)[[1]][1])
pos_trb_d_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[2],"*",fixed=T)[[1]][1])
pos_trb_j_gene <- sapply(strsplit(pos[,1],"@"),function(x) strsplit(x[3],"*",fixed=T)[[1]][1])
pos_df <- data.frame(TRB_cdr3=pos_trb_cdr3,TRB_v_gene=pos_trb_v_gene,
TRB_d_gene=pos_trb_d_gene,TRB_j_gene=pos_trb_j_gene,Clonality="na",
pmhc_code="EGFRpL858R",id="KITDFGRAK",binds=1,TCR_id=seq(nrow(pos)),
stringsAsFactors=F)
pos_df <- pos_df[grep("_", pos_df[,1], invert = T,fixed=T),]
pos_df <- pos_df[grep("*", pos_df[,1], invert = T,fixed=T),]
neg_trb_cdr3 <- sapply(strsplit(neg[,1],"@"),function(x) x[5])
neg_trb_v_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[1],"*",fixed=T)[[1]][1])
neg_trb_d_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[2],"*",fixed=T)[[1]][1])
neg_trb_j_gene <- sapply(strsplit(neg[,1],"@"),function(x) strsplit(x[3],"*",fixed=T)[[1]][1])
neg_df <- data.frame(TRB_cdr3=neg_trb_cdr3,TRB_v_gene=neg_trb_v_gene,
TRB_d_gene=neg_trb_d_gene,TRB_j_gene=neg_trb_j_gene,Clonality="na",
pmhc_code="EGFRpL858R",id="KITDFGRAK",binds=0,TCR_id=seq(nrow(neg))+nrow(pos),
stringsAsFactors=F)
neg_df <- neg_df[grep("_", neg_df[,1], invert = T,fixed=T),]
neg_df <- neg_df[grep("*", neg_df[,1], invert = T,fixed=T),]
tran_df <- rbind(pos_df,neg_df)
write.table(tran_df, file="EGFRpL858R.tran_df.mo3.csv",sep=",", quote=FALSE, row.names=F, col.names=T)

###########################################################
#building model by TCRAI framework (python version 3.6.15)#
###########################################################
#${TrainFile} is the file of training datasets
#${mo3File} is the valudation file from 3 month
#${mo1File} is the valudation file from 1 month
#${mo2Data} is the data from 2 month
#${mo2Label} is the label from 2 month
#${mo0Data} is the data from pretreatment
#${mo1Data} is the data from 1 month
#${mo3Data} is the data from 3 month

#training model
import os;import sys;import tempfile;import io;import numpy as np;import pandas as pd;
from sklearn.model_selection import train_test_split, ParameterGrid;
from sklearn.utils.class_weight import compute_class_weight;
from sklearn.model_selection import StratifiedShuffleSplit;
from sklearn import metrics;from sklearn.cluster import KMeans;import tensorflow as tf;
from tensorflow.keras.preprocessing.text import Tokenizer;
import tensorflow.keras as keras;
from PIL import Image;
import matplotlib.pyplot as plt;
import seaborn as sns;
sns.set(context='talk',palette='bright');
import umap;
from tcrai.modelling.processing import AASeqProcessor;
from tcrai.modelling import extractors, closers;
from tcrai.modelling.classification import SeqClassificationModelWithProcessor, VJXCRSeqModel;
from tcrai.plotting import ml_plots;
from tcrai.motif import dim_reduction,logo,motif_extraction

os.environ['PYTHONHASHSEED']=str(13);
np.random.seed(13);
tf.random.set_seed(13);
data_root = TrainFile;
df = pd.read_csv(data_root);df.id.unique(); df = df[df['id']=='KITDFGRAK'];df['labels'] = df['binds'];
df = df.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
max_len_cdr3 = 40;processor_b = AASeqProcessor(max_len_cdr3);
seq_processors = {
    'TRB_cdr3': processor_b
}
vj_cols = ['TRB_v_gene','TRB_j_gene']; vj_tokenizers = dict();
for col in vj_cols:
    vj_tokenizers[col] = Tokenizer(filters='',lower=False,oov_token='UNK')
    vj_tokenizers[col].fit_on_texts(df[col])
seq_encoder_params = {
    'embed_dim': 16,
    'filters': [64,128,256],
    'kernel_widths': [5,4,4],
    'dilations': [1,1,1],
    'strides': [1,3,3],
    'L2_conv': 0.01,
    'dropout_conv': 0.3,
};
seq_encoders = {
    'TRB_cdr3': extractors.conv_seq_extractor(seq_encoder_params,
                                              max_len_cdr3,
                                              seq_processors['TRB_cdr3'].vocab_size)
};
vj_encs = dict();
embed_dim = {
    'TRB_v_gene': 16,
    'TRB_j_gene': 8
};
for v in embed_dim.keys():
    hp_vj={
        'vj_width':len(vj_tokenizers[v].word_index) +1,
        'vj_embed': embed_dim[v],
        'dropout': 0.3
    }
    vj_encs[v] = extractors.vj_extractor(hp_vj,name=v+'_enc')
hp_pred = {
    'units': [],
    'dropout': 0.3,
    'L2_dense': 0.0,
    'init_bias': 0.0
};
predictor = closers.make_predictor(hp_pred);
model0 = VJXCRSeqModel(seq_encoders,
                       vj_encs,
                       predictor,
                       input_list=None 
                      )
model = SeqClassificationModelWithProcessor(model0, processors=seq_processors, extra_tokenizers=vj_tokenizers);
def df_to_input(df):
    """ convert a dataframe into a dictionary of inputs """
    cols = ['TRB_v_gene','TRB_j_gene','TRB_cdr3']
    x = { c:df[c].values for c in cols}
    return x
tr_df,test_df,y_tr,y_test = train_test_split(df,
                                            df['labels'],
                                           test_size=0.1,
                                           stratify=df['labels'],
                                           random_state=42)

train_df,vali_df,y_train,y_vali = train_test_split(tr_df,
                                            tr_df['labels'],
                                           test_size=0.2,
                                           stratify=tr_df['labels'],
                                           random_state=42)
train_input = df_to_input(train_df);
vali_input = df_to_input(vali_df);
test_input = df_to_input(test_df);
model.compile(
    optimizer = keras.optimizers.Adam(),
    loss = keras.losses.BinaryCrossentropy(),
    metrics = [keras.metrics.AUC(name='ROC')]
);
weights = compute_class_weight('balanced',classes=np.arange(2),y=np.squeeze(y_train));
print(weights);
early_stop = keras.callbacks.EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True )
history = model.fit(train_input, y_train,
          validation_data = (vali_input, y_vali),
            epochs = 250,
            callbacks = [early_stop],
            class_weight =  {i:weights[i] for i in range(len(weights))},
            batch_size = 256,
            verbose = 0
);
model.add_evaluator(metrics.roc_auc_score,'roc');
model.evaluate(test_input,y_test)

#validation of datasets from 3 month
data_root3 = mo3File
df3 = pd.read_csv(data_root3);
df3.id.unique(); 
df3 = df3[df3['id']=='KITDFGRAK'];
df3['labels'] = df3['binds'];
df3 = df3.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])

test_df3= df3[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
y_test3 = df3['labels']										   
test_input3 = df_to_input(test_df3);
model.evaluate(test_input3,y_test3)

#validation of datasets from 1 month
data_root1 = mo1File
df1 = pd.read_csv(data_root1);
df1.id.unique(); 
df1 = df1[df1['id']=='KITDFGRAK'];
df1['labels'] = df1['binds'];
df1 = df1.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
test_df1= df1[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
y_test1 = df1['labels']
test_input1 = df_to_input(test_df1);
model.evaluate(test_input1,y_test1)
preds = model.run(test_input)

#plot the ROC curve
ml_plots.plot_roc_binomial(y_test,preds,save_path='mo2.roc.pdf')
preds3 = model.run(test_input3)
ml_plots.plot_roc_binomial(y_test3,preds3,save_path='mo3.roc.pdf')
preds1 = model.run(test_input1)
ml_plots.plot_roc_binomial(y_test1,preds1,save_path='mo1.roc.pdf')
predsT = model.run(train_input)
ml_plots.plot_roc_binomial(y_train,predsT,save_path='mo2.train.roc.pdf')
predsV = model.run(vali_input)
ml_plots.plot_roc_binomial(y_vali,predsV,save_path='mo2.vali.roc.pdf')
train_df.to_csv(mo2Data, index=False, header=True )
y_train.to_csv(mo2Label, index=False, header=True )
np.savetxt('mo2.train.predict.txt', predsT, fmt="%.4f", delimiter=',')

#pretreatment
data_all0 = mo0Data
dfa0 = pd.read_csv(data_all0);
dfa0.id.unique(); 
dfa0 = dfa0[dfa0['id']=='KITDFGRAK'];
dfa0['labels'] = dfa0['binds'];
dfa0 = dfa0.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
test_dfa0= dfa0[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
y_test_dfa0 = dfa0['labels']
test_input_dfa0 = df_to_input(test_dfa0);
preds_dfa0 = model.run(test_input_dfa0)
np.savetxt('mo0.all.predict.txt', preds_dfa0, fmt="%.4f", delimiter=',')
dfa0.to_csv('mo0.all.process.df.txt', index=False, header=True )

#1 month
data_all1 = mo1Data
dfa1 = pd.read_csv(data_all1);
dfa1.id.unique(); 
dfa1 = dfa1[dfa1['id']=='KITDFGRAK'];
dfa1['labels'] = dfa1['binds'];
dfa1 = dfa1.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
test_dfa1= dfa1[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
y_test_dfa1 = dfa1['labels']
test_input_dfa1 = df_to_input(test_dfa1);
preds_dfa1 = model.run(test_input_dfa1)
np.savetxt('mo1.all.predict.txt', preds_dfa1, fmt="%.4f", delimiter=',')
dfa1.to_csv('mo1.all.process.df.txt', index=False, header=True )

#2 month
data_all2 = mo2Data
dfa2 = pd.read_csv(data_all2);
dfa2.id.unique(); 
dfa2 = dfa2[dfa2['id']=='KITDFGRAK'];
dfa2['labels'] = dfa2['binds'];
dfa2 = dfa2.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
test_dfa2= dfa2[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
y_test_dfa2 = dfa2['labels']
test_input_dfa2 = df_to_input(test_dfa2);
preds_dfa2 = model.run(test_input_dfa2)
np.savetxt('mo2.all.predict.txt', preds_dfa2, fmt="%.4f", delimiter=',')
dfa2.to_csv('mo2.all.process.df.txt', index=False, header=True )

#3 month
data_all3 = mo3Data
dfa3 = pd.read_csv(data_all3);
dfa3.id.unique(); 
dfa3 = dfa3[dfa3['id']=='KITDFGRAK'];
dfa3['labels'] = dfa3['binds'];
dfa3 = dfa3.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
test_dfa3= dfa3[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
y_test_dfa3 = dfa3['labels']
test_input_dfa3 = df_to_input(test_dfa3);
preds_dfa3 = model.run(test_input_dfa3)
np.savetxt('mo3.all.predict.txt', preds_dfa3, fmt="%.4f", delimiter=',')
dfa3.to_csv('mo3.all.process.df.txt', index=False, header=True )

###########################################
#apply model to CD8+ T cells from scRNAseq#
###########################################
#${smp} is the sample ID
#${workDir} is the working directory
#${ImmuneFile} is the seurat obeject for all the immune cell
#${TRUSTdir} is the directory of TRUST result
#${TRUSTres} is the result of TRUST, which formated ad input for model
#${ClinicalFile} is the file of clinical information
#${immuneFile} is the seurat object of all the immune cells

#TRUST analysis
run-trust4 -b ${smp}.bam -f hg38_bcrtcr.fa --ref human_IMGT+C.fa --barcode CB -o LSF_L -t 16

#make the input for model
#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(gplots)

#Load all the immune cell
immu.integrated <- readRDS(file = ImmuneFile)

#Load the results of TRUST
setwd(TRUSTdir)
fl <- dir()
fl <- grep("_barcode_report.tsv",fl,value=T)
fl_list <- sapply(fl,function(x){
	fl_tmp <- read.table(file =x,sep="\t",stringsAsFactors=F, header=F)
	fl_tmp <- fl_tmp[(fl_tmp[,2]=="abT")&(fl_tmp[,3]!="*"),]
	fl_tmp_df <- data.frame(V1=fl_tmp[,1],
	V2=sapply(strsplit(fl_tmp[,3],","),function(y) y[1]),
	V3=sapply(strsplit(fl_tmp[,3],","),function(y) y[3]),
	V4=sapply(strsplit(fl_tmp[,3],","),function(y) y[6]),stringsAsFactors=F)
	fl_tmp_df[,2] <- sapply(strsplit(fl_tmp_df[,2],"*",fixed=T),function(y) y[1])
	fl_tmp_df[,3] <- sapply(strsplit(fl_tmp_df[,3],"*",fixed=T),function(y) y[1])
	smp <- sub("_barcode_report.tsv", "", x)
	print(smp)
	anno <- subset(immu.integrated,
	cells=rownames(immu.integrated@meta.data)[(as.character(immu.integrated@meta.data[,"orig.ident"])==smp)&(as.character(immu.integrated@meta.data[,"majorLabel"])=="CD8T")])
	rownames(anno@meta.data) <- sapply(strsplit(rownames(anno@meta.data),"_"),function(y) y[1])
	fl_tmp_df <- fl_tmp_df[fl_tmp_df[,1]%in%rownames(anno@meta.data),]
	label <- paste(fl_tmp_df[,4],fl_tmp_df[,2],fl_tmp_df[,3],sep="_")
	fl_tmp_label_df <- cbind(fl_tmp_df,label)
	fl_tmp_label_df <- fl_tmp_label_df[!duplicated(fl_tmp_label_df[,"label"]),]
	if(nrow(fl_tmp_label_df)>0)
	{
		nclone <- sapply(fl_tmp_label_df[,"label"],function(y) sum(label==y))
		freq <- nclone/sum(nclone)
		res_df <- data.frame(TRB_cdr3=fl_tmp_label_df[,"V4"],
		TRB_v_gene=fl_tmp_label_df[,"V2"],TRB_j_gene=fl_tmp_label_df[,"V3"],
		Clonality=freq,pmhc_code="EGFRpL858R",id="KITDFGRAK",Sample=smp,
		stringsAsFactors=F)
		rownames(res_df) <- NULL
		return(res_df)
	}
})
fl_df <- do.call(rbind,fl_list)
rownames(fl_df) <- NULL
TCR_id <- seq(nrow(fl_df))
fl_df <- cbind(fl_df,TCR_id)

#Set working directory
setwd(workDir)
write.table(fl_df, file=TRUSTres, sep=",",quote=FALSE, row.names=F, col.names=T )  

#apply model to CD8+ T cells from scRNAseq
data_our = TRUSTres
dfaour = pd.read_csv(data_our);
dfaour.id.unique(); 
dfaour = dfaour[dfaour['id']=='KITDFGRAK'];
dfaour = dfaour.drop_duplicates(subset=['TRB_v_gene',
                                'TRB_j_gene',
                                'TRB_cdr3'
                               ])
test_dfaour= dfaour[['TRB_v_gene','TRB_j_gene','TRB_cdr3']]
test_input_dfaour = df_to_input(test_dfaour);
preds_dfaour = model.run(test_input_dfaour)
np.savetxt('dfaour.predict.txt', preds_dfaour, fmt="%.4f", delimiter=',')
dfaour.to_csv('dfaour.process.df.txt', index=False, header=True )

#identification of differetially expression gene of neoantigen-recognizing T cells 
#among clinical stages

#Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(gplots)

#Load the clinical information and the model result
clinical <- read.table(file =ClinicalFile,sep="\t",quote = "",stringsAsFactors=F, header=T)
clinical[clinical[,"class"]%in%c("LUSC","Dermoid","Sarcomatoid"),"class"] <- "NonAden"
clinical_label <- paste(clinical[,5],clinical[,4],sep="_")
TCRAIinput <- read.table(file = "ourCD8T.csv",sep=",",quote = "",stringsAsFactors=F,header=T)
TCRAIres <- read.table(file = "dfaour.process.df.txt",sep=",",quote = "",stringsAsFactors=F,header=T)
TCRAIpredict <- read.table(file = "dfaour.predict.txt",sep=",",quote = "",stringsAsFactors=F,header=F)
TCRAIres_predict <- cbind(TCRAIres,TCRAIpredict)
inputlabel <- paste(TCRAIinput[,"TRB_cdr3"],TCRAIinput[,"TRB_v_gene"],TCRAIinput[,"TRB_j_gene"],sep="_")
prelabel <- paste(TCRAIres_predict[,"TRB_cdr3"],TCRAIres_predict[,"TRB_v_gene"],
TCRAIres_predict[,"TRB_j_gene"],sep="_")
inputPredict <- TCRAIres_predict[match(inputlabel,prelabel),"V1"]
TCRAIinput_predict <- cbind(TCRAIinput,inputPredict)
immu.integrated <- readRDS(file = immuneFile)

#Load the results of TRUST
setwd(TRUSTdir)
fl <- dir()
fl <- grep("_barcode_report.tsv",fl,value=T)
fl_list <- sapply(fl,function(x){
	fl_tmp <- read.table(file =x,sep="\t",stringsAsFactors=F, header=F)
	fl_tmp <- fl_tmp[(fl_tmp[,2]=="abT")&(fl_tmp[,3]!="*"),]
	fl_tmp_df <- data.frame(V1=fl_tmp[,1],
	V2=sapply(strsplit(fl_tmp[,3],","),function(y) y[1]),
	V3=sapply(strsplit(fl_tmp[,3],","),function(y) y[3]),
	V4=sapply(strsplit(fl_tmp[,3],","),function(y) y[6]),stringsAsFactors=F)
	fl_tmp_df[,2] <- sapply(strsplit(fl_tmp_df[,2],"*",fixed=T),function(y) y[1])
	fl_tmp_df[,3] <- sapply(strsplit(fl_tmp_df[,3],"*",fixed=T),function(y) y[1])
	smp <- sub("_barcode_report.tsv", "", x)
	print(smp)
	anno <- subset(immu.integrated,
	cells=rownames(immu.integrated@meta.data)[(as.character(immu.integrated@meta.data[,"orig.ident"])==smp)&(as.character(immu.integrated@meta.data[,"majorLabel"])=="CD8T")])
	smp_index <- strsplit(rownames(anno@meta.data)[1],"_")[[1]][2]
	fl_tmp_df[,1] <- paste(fl_tmp_df[,1],smp_index,sep="_")
	fl_tmp_df <- fl_tmp_df[fl_tmp_df[,1]%in%rownames(anno@meta.data),]

	if(nrow(fl_tmp_df)>0)
	{
		label <- paste(fl_tmp_df[,4],fl_tmp_df[,2],fl_tmp_df[,3],sep="_")
		fl_tmp_label_df <- cbind(fl_tmp_df,label)
		cell_state <- anno@meta.data[fl_tmp_df[,1],"minorLabel"]
		fl_tmp_label_df <- cbind(fl_tmp_label_df,cell_state)
		res_df <- data.frame(TRB_cdr3=fl_tmp_label_df[,"V4"],
		TRB_v_gene=fl_tmp_label_df[,"V2"],TRB_j_gene=fl_tmp_label_df[,"V3"],
		pmhc_code="EGFRpL858R",id="KITDFGRAK",Sample=smp,cell_state=cell_state,
		barcode=fl_tmp_label_df[,"V1"],
		stringsAsFactors=F)
		rownames(res_df) <- NULL
		return(res_df)
	}
})

fl_df <- do.call(rbind,fl_list)
rownames(fl_df) <- NULL
fl_df_label <- paste(fl_df[,1],fl_df[,2],fl_df[,3],sep="_")
predict_v1 <- TCRAIinput_predict[match(fl_df_label,inputlabel),"inputPredict"]
fl_df <- cbind(fl_df,predict_v1)
fl_df_pos <- fl_df[fl_df[,"predict_v1"]>0.81595,]
EGFR_smp_list <- c("WZH","HKJ","ZQZ","XBC","ZL","LTJ","WHZ","YFQ","WLB","CXM","BYQ",
"LYQ","WXY","CSM")
#samples with EGFR L858R neoantigen
fl_df_pos_bind <- fl_df_pos[fl_df_pos[,"Sample"]%in%EGFR_smp_list,]
group <- clinical[match(fl_df_pos_bind[,"Sample"],clinical[,"rawDataID"]),"class"]
fl_df_pos_bind <- cbind(fl_df_pos_bind,group)
TCR_bind <- subset(immu.integrated,cells=fl_df_pos_bind[,"barcode"])
TCR_bind$"clinical" <-  fl_df_pos_bind[match(rownames(TCR_bind@meta.data),
fl_df_pos_bind[,"barcode"]),"group"]
DefaultAssay(TCR_bind)<-"RNA"

#find markers
CD8T.markers_early<-FindMarkers(TCR_bind,ident.1=c("HiDenGGO","pGGO"),
ident.2=c("s25GGO","s50GGO","s75GGO","Solid1","SolidN"),
group.by = "clinical")
CD8T.markers_early1<-FindMarkers(TCR_bind,ident.1=c("HiDenGGO","pGGO","s25GGO","s50GGO"),
ident.2=c("s75GGO","Solid1","SolidN"),
group.by = "clinical")
CD8T.markers_early2<-FindMarkers(TCR_bind,ident.1=c("HiDenGGO","pGGO","s25GGO"),
ident.2=c("s75GGO","Solid1","SolidN"),
group.by = "clinical")
CD8T.markers_early3<-FindMarkers(TCR_bind,ident.1=c("pGGO","s25GGO","s50GGO"),
ident.2=c("s75GGO","Solid1","SolidN"),
group.by = "clinical")
CD8T.markers_early4<-FindMarkers(TCR_bind,ident.1=c("HiDenGGO","pGGO"),
ident.2=c("s25GGO","s50GGO"),
group.by = "clinical")
#s25GGO,s50GGO,s75GGO and HiDenGGO are also called 
#GGO25, GGO50,GGO75 and dGGO in the manuscripts