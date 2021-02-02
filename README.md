# CODEX-KNN

CODEX-KNN is the new pipeline for calling CNVs based on the CODEX tool, where the normalization process is conducted separately for each of the sequencing regions. Independent normalization allows to significantly improve the results of rare CNVs detection. For example, for investigated dataset we reduced the number of false positives calls from over 15,000 to around 5,000 while maintaining a constant number of true positive calls equal to about 150 CNVs. However, independent normalization of each sequencing region is a computationally expensive process, therefore our pipeline is customized and can be easily run in the cloud computing environment, on the computer cluster or on the single CPU server. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

You need to install the below software:

```
Docker
```

## Running the tests

The CODEX-KNN tool is composed of several steps, below we've got step-by-step instructions how to run the application. Before starting the tests below, please copy the contents of the dataset and scripts directories to the /tmp directory.
```
cd CODEX-KNN
cp dataset/* /tmp
cp scripts/* /tmp
cd /tmp
```

### Quality control

The first step of the CODEX-KNN pipeline is quality control. In this process all sequencing regions (I) with GC content below 20% or greater than 80%, (II) with mappability factor below 0.9, (III) with median read depth across all samples below 20 or greater than 4000 and (IV) shorter than 20 bp or longer than 2000 bp are removed. 

```
docker run --rm -v /tmp:/tmp -w /tmp biodatageeks/cnv-opt-target-qc Rscript -e "library('TARGET.QC');run_TARGET.QC(0.9,20,4000,20,2000,20,80,'/tmp/raw_cov.csv','/tmp/raw_cov_qc.csv','/tmp/exons_coordinates.bed','/tmp/exons_coordinates_qc.bed')"
```

### Normalization

After quality control the independent normalization process is applied. In our example we've got 100 sequencing regions, here, we parallelize and divide the calculations into four threads, in the normalization process we used for each sequencing region 50 most correlated background-modeling sequencing regions. After normalization we join the results and add the header to the output file.

```
for V in {1..25}
do
  docker run --rm -v /tmp:/tmp -v /tmp:/tmp -w /tmp biodatageeks/cnv-opt-codexcov Rscript -e "library('CODEXCOV');source('/tmp/run_CODEXCOV_knn_by_regions.R');run_CODEXCOV(1, 3, 200, 50, $V, 'raw_cov_qc.csv', 'exons_coordinates_qc.bed', '/tmp/norm_cov_qc_$V.csv')"
done &
for V in {26..50}
do
  docker run --rm -v /tmp:/tmp -v /tmp:/tmp -w /tmp biodatageeks/cnv-opt-codexcov Rscript -e "library('CODEXCOV');source('/tmp/run_CODEXCOV_knn_by_regions.R');run_CODEXCOV(1, 3, 200, 50, $V, 'raw_cov_qc.csv', 'exons_coordinates_qc.bed', '/tmp/norm_cov_qc_$V.csv')"
done &
for V in {51..75}
do
  docker run --rm -v /tmp:/tmp -v /tmp:/tmp -w /tmp biodatageeks/cnv-opt-codexcov Rscript -e "library('CODEXCOV');source('/tmp/run_CODEXCOV_knn_by_regions.R');run_CODEXCOV(1, 3, 200, 50, $V, 'raw_cov_qc.csv', 'exons_coordinates_qc.bed', '/tmp/norm_cov_qc_$V.csv')"
done &
for V in {76..100}
do
  docker run --rm -v /tmp:/tmp -v /tmp:/tmp -w /tmp biodatageeks/cnv-opt-codexcov Rscript -e "library('CODEXCOV');source('/tmp/run_CODEXCOV_knn_by_regions.R');run_CODEXCOV(1, 3, 200, 50, $V, 'raw_cov_qc.csv', 'exons_coordinates_qc.bed', '/tmp/norm_cov_qc_$V.csv')"
done &

ls -v /tmp/norm_cov_qc_{1..100}.csv | xargs cat > /tmp/norm_cov_qc.csv

echo -e "NA06984,NA06985,NA06986,NA06989,NA06994,NA07000,NA07037,NA07048,NA07051,NA07056,NA07347,NA07357,NA10847,NA10851,NA11829,NA11830,NA11831,NA11832,NA11840,NA11843,NA11881,NA11892,NA11893,NA11894,NA11918,NA11919,NA11920,NA11930,NA11931,NA11932,NA11933,NA11992,NA11994,NA11995,NA12003,NA12004,NA12005,NA12006,NA12043,NA12044,NA12045,NA12046,NA12058,NA12144,NA12154,NA12155,NA12156,NA12234,NA12249,NA12272,NA12273,NA12275,NA12282,NA12283,NA12286,NA12287,NA12340,NA12341,NA12342,NA12347,NA12348,NA12383,NA12399,NA12400,NA12413,NA12414,NA12489,NA12546,NA12716,NA12717,NA12718,NA12748,NA12749,NA12750,NA12751,NA12760,NA12761,NA12762,NA12763,NA12775,NA12776,NA12777,NA12778,NA12812,NA12813,NA12814,NA12815,NA12827,NA12828,NA12829,NA12830,NA12842,NA12843,NA12872,NA12873,NA12874,NA12878,NA12889,NA12890,NA18486,NA18488,NA18489,NA18498,NA18499,NA18501,NA18502,NA18504,NA18505,NA18507,NA18508,NA18510,NA18511,NA18516,NA18517,NA18519,NA18520,NA18522,NA18523,NA18525,NA18526,NA18528,NA18530,NA18531,NA18532,NA18533,NA18534,NA18535,NA18536,NA18537,NA18538,NA18539,NA18541,NA18542,NA18543,NA18544,NA18545,NA18546,NA18547,NA18548,NA18549,NA18550,NA18552,NA18553,NA18555,NA18557,NA18558,NA18559,NA18560,NA18561,NA18562,NA18563,NA18564,NA18565,NA18566,NA18567,NA18570,NA18571,NA18572,NA18573,NA18574,NA18577,NA18579,NA18582,NA18591,NA18592,NA18593,NA18595,NA18596,NA18597,NA18599,NA18602,NA18603,NA18605,NA18606,NA18608,NA18609,NA18610,NA18611,NA18612,NA18613,NA18614,NA18615,NA18616,NA18617,NA18618,NA18619,NA18620,NA18621,NA18622,NA18623,NA18624,NA18625,NA18626,NA18627,NA18628,NA18629,NA18630,NA18631,NA18632,NA18633,NA18634,NA18635,NA18636,NA18637,NA18638,NA18639,NA18640,NA18641,NA18642,NA18643,NA18644,NA18645,NA18646,NA18647,NA18648,NA18740,NA18745,NA18747,NA18748,NA18749,NA18757,NA18853,NA18856,NA18858,NA18861,NA18864,NA18865,NA18867,NA18868,NA18870,NA18871,NA18873,NA18874,NA18876,NA18877,NA18878,NA18879,NA18881,NA18907,NA18908,NA18909,NA18910,NA18912,NA18915,NA18916,NA18917,NA18923,NA18924,NA18933,NA18934,NA18939,NA18940,NA18941,NA18942,NA18943,NA18944,NA18945,NA18946,NA18947,NA18948,NA18949,NA18950,NA18951,NA18952,NA18953,NA18954,NA18956,NA18957,NA18959,NA18960,NA18961,NA18962,NA18963,NA18964,NA18965,NA18966,NA18967,NA18968,NA18969,NA18970,NA18971,NA18972,NA18973,NA18974,NA18975,NA18976,NA18977,NA18978,NA18979,NA18980,NA18981,NA18982,NA18983,NA18984,NA18985,NA18986,NA18987,NA18988,NA18989,NA18990,NA18991,NA18992,NA18993,NA18994,NA18995,NA18997,NA18998,NA18999,NA19000,NA19001,NA19002,NA19003,NA19004,NA19005,NA19006,NA19007,NA19009,NA19010,NA19011,NA19012,NA19017,NA19019,NA19020,NA19023,NA19024,NA19025,NA19026,NA19027,NA19028,NA19030,NA19031,NA19035,NA19036,NA19037,NA19038,NA19041,NA19042,NA19043,NA19054,NA19055,NA19056,NA19057,NA19058,NA19059,NA19060,NA19062,NA19063,NA19064,NA19065,NA19066,NA19067,NA19068,NA19070,NA19072,NA19074,NA19075,NA19076,NA19077,NA19078,NA19079,NA19080,NA19081,NA19082,NA19083,NA19084,NA19085,NA19086,NA19087,NA19088,NA19089,NA19090,NA19091,NA19092,NA19093,NA19095,NA19096,NA19098,NA19099,NA19102,NA19107,NA19108,NA19113,NA19114,NA19116,NA19117,NA19118,NA19119,NA19121,NA19129,NA19130,NA19131,NA19137,NA19138,NA19141,NA19143,NA19144,NA19146,NA19147,NA19149,NA19152,NA19153,NA19159,NA19160,NA19171,NA19172,NA19175,NA19184,NA19185,NA19189,NA19190,NA19197,NA19198,NA19200,NA19201,NA19204,NA19206,NA19207,NA19209,NA19210,NA19213,NA19214,NA19222,NA19223,NA19225,NA19235,NA19236,NA19238,NA19239,NA19240,NA19247,NA19248,NA19256,NA19257,NA19307,NA19308,NA19309,NA19310,NA19311,NA19312,NA19313,NA19314,NA19315,NA19316,NA19317,NA19318,NA19319,NA19320,NA19321,NA19323,NA19324,NA19327,NA19328,NA19331,NA19332,NA19334,NA19338,NA19346,NA19347,NA19350,NA19351,NA19355,NA19360,NA19372,NA19374,NA19375,NA19376,NA19377,NA19378,NA19379,NA19380,NA19383,NA19384,NA19385,NA19390,NA19391,NA19393,NA19394,NA19395,NA19397,NA19399,NA19401,NA19403,NA19404,NA19428,NA19429,NA19430,NA19431,NA19434,NA19435,NA19436,NA19437,NA19438,NA19439,NA19440,NA19443,NA19445,NA19446,NA19448,NA19449,NA19451,NA19452,NA19454,NA19455,NA19456,NA19457,NA19461,NA19462,NA19463,NA19466,NA19467,NA19468,NA19471,NA19472,NA19473,NA19474,NA19475,NA19625,NA19648,NA19649,NA19651,NA19652,NA19654,NA19655,NA19657,NA19658,NA19660,NA19661,NA19663,NA19664,NA19669,NA19670,NA19675,NA19676,NA19678,NA19679,NA19681,NA19682,NA19684,NA19685,NA19700,NA19701,NA19703,NA19704,NA19707,NA19711,NA19712,NA19713,NA19716,NA19717,NA19719,NA19720,NA19722,NA19723,NA19725,NA19726,NA19728,NA19729,NA19731,NA19732,NA19734,NA19735,NA19740,NA19741,NA19746,NA19747,NA19749,NA19750,NA19752,NA19755,NA19756,NA19758,NA19759,NA19761,NA19762,NA19764,NA19770,NA19771,NA19773,NA19774,NA19776,NA19777,NA19779,NA19780,NA19782,NA19783,NA19785,NA19786,NA19788,NA19789,NA19792,NA19794,NA19795,NA19818,NA19819,NA19834,NA19835,NA19900,NA19901,NA19904,NA19908,NA19909,NA19913,NA19914,NA19916,NA19917,NA19920,NA19921,NA19922,NA19923,NA19982,NA19984,NA19985,NA20126,NA20127,NA20274,NA20276,NA20278,NA20281,NA20282,NA20287,NA20289,NA20291,NA20294,NA20296,NA20298,NA20299,NA20314,NA20317,NA20318,NA20320,NA20321,NA20322,NA20332,NA20334,NA20336,NA20339,NA20340,NA20341,NA20342,NA20344,NA20346,NA20348,NA20351,NA20355,NA20356,NA20357,NA20359,NA20362,NA20412,NA20502,NA20503,NA20504,NA20505,NA20506,NA20507,NA20508,NA20509,NA20510,NA20511,NA20512,NA20513,NA20514,NA20515,NA20516,NA20517,NA20518,NA20519,NA20520,NA20521,NA20522,NA20524,NA20525,NA20526,NA20527,NA20528,NA20529,NA20530,NA20531,NA20532,NA20533,NA20534,NA20535,NA20536,NA20538,NA20539,NA20540,NA20541,NA20542,NA20543,NA20544,NA20581,NA20582,NA20585,NA20586,NA20587,NA20588,NA20589,NA20752,NA20753,NA20754,NA20755,NA20756,NA20757,NA20758,NA20759,NA20760,NA20761,NA20762,NA20763,NA20764,NA20765,NA20766,NA20767,NA20768,NA20769,NA20770,NA20771,NA20772,NA20773,NA20774,NA20775,NA20778,NA20783,NA20785,NA20786,NA20787,NA20790,NA20792,NA20795,NA20796,NA20797,NA20798,NA20799,NA20802,NA20803,NA20804,NA20805,NA20806,NA20807,NA20808,NA20809,NA20810,NA20811,NA20812,NA20813,NA20814,NA20815,NA20818,NA20819,NA20821,NA20822,NA20826,NA20827,NA20828,NA20832,NA20845,NA20846,NA20847,NA20849,NA20850,NA20851,NA20852,NA20853,NA20854,NA20856,NA20858,NA20859,NA20861,NA20862,NA20863,NA20864,NA20866,NA20867,NA20868,NA20869,NA20870,NA20871,NA20872,NA20874,NA20875,NA20876,NA20877,NA20878,NA20881,NA20882,NA20884,NA20885,NA20886,NA20887,NA20888,NA20889,NA20890,NA20891,NA20892,NA20893,NA20894,NA20895,NA20896,NA20897,NA20898,NA20899,NA20900,NA20901,NA20902,NA20903,NA20904,NA20905,NA20906,NA20908,NA20910,NA20911,NA21086,NA21087,NA21088,NA21089,NA21090,NA21091,NA21092,NA21093,NA21094,NA21095,NA21097,NA21098,NA21099,NA21100,NA21101,NA21102,NA21103,NA21104,NA21105,NA21106,NA21107,NA21108,NA21109,NA21110,NA21111,NA21112,NA21113,NA21114,NA21115,NA21116,NA21117,NA21118,NA21119,NA21120,NA21122,NA21123,NA21124,NA21125,NA21126,NA21127,NA21128,NA21129,NA21130,NA21133,NA21135,NA21137,NA21141,NA21142,NA21143,NA21144" | cat - /tmp/norm_cov_qc.csv > /tmp/tmp.csv
mv /tmp/tmp.csv /tmp/norm_cov_qc.csv
```

### CNV calling

Finally, call the resultant set of CNVs with command:

```
docker run --rm -it -v /tmp:/tmp -w /tmp biodatageeks/cnv-opt-codexcov Rscript -e "library('CODEXCOV');source('/tmp/run_CODEXCOV_from_Yhat.R');run_CODEXCOV(1, 3, 200, '/tmp/raw_cov_qc.csv', '/tmp/norm_cov_qc.csv', '/tmp/exons_coordinates_qc.bed', '/tmp/reference_sample_set_kmeans_1.csv', '/tmp/calls.csv')"
```
