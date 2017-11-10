# SCRUB

Single Cell RNA batch effect correcting and biology variance detection

It's a droplet based single cell RNA batch effect correcting protocal....

One of the critical parts is the between batches variant genes detected, which is based on 2 statistical test and one similarity metric. 

It's basically a combination of the model free test (Mann-Whiteney test), model based test (KS test, Specifically, emperical distribution based) and Battachayya distance (the similarity metric of overlap of pdf).

