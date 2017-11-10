# SCRUB

Single Cell RNA batch effect correcting and biology variance detection method: SCRUB (Single Cell RNA Unified molecular identifier Batch effect correcting)

It's a droplet based single cell RNA batch effect correcting protocal....

One of the critical parts is the between batches variant genes detection, which is based on 2 statistical tests and one similarity metric. 

It's basically a combination of the model free test (Mann-Whiteney test), model based test (KS test, Specifically, emperical distribution based) and Battachayya distance (the similarity metric of overlap of pdf).

Also, a non-trivial expansion to multiple batches biology variance detection and batch effect correcting.
