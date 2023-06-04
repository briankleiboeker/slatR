# slatR
Shotgun Lipidomic Analysis Toolkit (an R package)

An R package version of [this webapp](https://github.com/briankleiboeker/SLAT) to enable high-througput, reprodicuble analysis and annotation of shotgun lipidomic data with some added functionality for making comparisons across biological conditions.  

### Install and load slatR
```{r}
install.packages("devtools")
devtools::install_github('briankleiboeker/slatR')

library(slatR)
```
### Read in data and assign elemental composition to m/z values
```
df <- read.csv('fulldata_clean.csv')

assignment_df <- assign_comp_from_mz(df$mz,
                               elements = c("C","C","H","O","N","P","Na","Cl"),
                               isotope = c(12,13,1,16,14,31,23,35),
                               mins.list = c(10,0,20,2,0,0,0,0),
                               maxs.list = c(60,5,160,16,2,2,0,1),
                               ionmode = "neg.ion",
                               error.ppm = 3,
                               rdb.rule = "nonint")
                               
df <- cbind(df,assignment_df[,-1])
```

### Prepare data: filter out lipids which are not detected in at least 1/2 of replicates in at least one biological condition
```{r}
library(tidyr)
library(dplyr)

df2 <- df %>% 
  mutate(condition = stringr::str_sub(fullsample,start = 1, end = 2)) %>% 
  group_by(composition,condition) %>% 
  mutate(num = n()) %>% 
  ungroup() %>% 
  group_by(composition) %>% 
  filter(max(num)>3) %>% 
  ungroup()
```
### Prepare data: Convert dataframe to a matrix in the 'wide' format (i.e. with rows corresponding to lipids/species, columns corresponding to samples/replicates)
```{r}
wide.matrix <- df2 %>% dplyr::select(c("rel.intensity","fullsample","composition")) %>% 
  filter(rel.intensity != 100) %>% 
  dplyr::rename("sample" = "fullsample",
                "comp" = "composition") %>% 
  pivot_wider(names_from = sample,
              values_from = rel.intensity,
              values_fn = mean)

wide.matrix <- as.matrix(wide.matrix)
rownames(wide.matrix)<-wide.matrix[,1]
wide.matrix<-wide.matrix[,-1]
```
### Normalize data 
```{r}
norm.mat <- normalize_sl(wide.matrix)
```
### Impute missing values 
```{r}
full.mat <- impute_sl(norm.mat)
```
### Test for differentially abundant lipids
```{r}
sc<-grep("sc",colnames(full.mat),value = T)
ko<-grep("ko",colnames(full.mat),value = T)

result.df <- test_for_da_sl(full.mat,sc,ko)
```
### Assign general lipid structures
```{r}
result.df <- result.df %>% tibble::rownames_to_column() %>%
             dplyr::rename("composition" = "rowname")
result.df$gen.structures <- assign_structures(comp = result.df$composition,
                                              "neutral",
                                              domain = "euk",
                                              max.dbl.bnds = 8)
```
### Guess exact lipid structures from exact structure database
```{r}
es.df <- guess_exact_structures(result.df$gen.structures)
results <- cbind(result.df,es.df)
```
### Visualize results!
![volcanoplot](https://github.com/briankleiboeker/slatR/assets/59810795/725b6df6-081d-459e-8387-173c7307690f)
