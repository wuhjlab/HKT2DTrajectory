## This code is used to search for all linear trajectories that meet the required number of people.

## df_directivity is the output from directivity tests.

subset_vector <- function(n,vec) {
  if (n > length(vec)) {
    return(vec)
  } else {
    return(vec[1:n])
  }
}

grep_two<-function(x){
  x1<-str_split_1(x,"\\_")
  x2<-lag(x1)[-1]
  x1<-x1[-1]
  y<-data.frame(node1=x2,node2=x1)%>%mutate(path=x)
  return(y)
}



Npop <- length(unique(dat$id))
res<-pblapply(n_last_edge:(n_last_edge+50),function(i){

  pvalue <- 1-sum(dbinom(n_last_edge:Npop, size=Npop, prob=i/Npop))
  res=data.frame(N=i,p=pvalue)
  return(res)
})
n_min<-bind_rows(res)%>%filter(p<0.001)
n_min<-n_min[1,1]


node<-df_directivity%>%filter(directivity)%>%distinct()

dat_dia<-dat%>%select(id,icd_main,date)%>%
  filter(icd_main%in%unique(c(node$exposure,node$outcome)))%>%
  group_by(id,icd_main)%>%filter(date==min(date))%>%distinct()%>%ungroup()%>%
  arrange(id,date)



index<-unique(dat_dia$id)

dat_path<-mclapply(index,function(index){

  data<-dat_dia%>%filter(id==index)
  res<-get_path(data)
  return(res)
},mc.cores=thread)

dat_path<- as.data.frame(rbindlist(dat_path))



dat_path_num<-dat_path%>%filter(grepl("\\_",path))%>%select(-path)%>%distinct()%>%group_by(sub_path)%>%
  summarise(freq=n())%>%filter(freq>=n_min)
dat_path_num$path_length<-str_count(dat_path_num$sub_path,fixed("_"))+1

g<-graph_from_data_frame(node%>%select(exposure,outcome,RR),directed = T)

get_paths_for_vertex <- function(i) {
  paths_for_i <- list()
  for (j in 1:vcount(g)) {
    if (i != j) {
      paths <- all_simple_paths(g, from=V(g)[i], to=V(g)[j])
      paths_for_i <- c(paths_for_i, paths)
    }
  }
  return(paths_for_i)
}

all_paths<- mclapply(1:vcount(g), get_paths_for_vertex, mc.cores=thread)
all_paths<- unlist(all_paths,
                   recursive = FALSE)

all_paths_name<-mclapply(all_paths,function(x)paste(V(g)[x]$name, collapse = "_"),mc.cores = thread)
all_paths_name<-data.frame(sub_path=unlist(all_paths_name))

dat_path_num_final<-merge(all_paths_name,dat_path_num,by='sub_path')

dat_sub_path_name<-rbindlist(mclapply(dat_path_num_final$sub_path,grep_subgraph,mc.cores=thread))
res_final<-dat_sub_path_name%>%left_join(dat_path_num_final,by='sub_path',
                                         relationship = "many-to-many")%>%
  rename('patients_number'=freq)

dat_two<-rbindlist(mclapply(dat_path_num_final$sub_path,grep_two,mc.cores = thread))
dat_two<-dat_two%>%left_join(node%>%
                               rename('node1'=exposure,'node2'=outcome),
                             by=c('node1', 'node2'),relationship = "many-to-many")

res_final<-res_final%>%left_join(dat_two%>%rename('sub_path'=path),by='sub_path')
res_final$path_length<-str_count(res_final$path,fixed("_"))+1


