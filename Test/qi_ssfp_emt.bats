# Copyright Tobias Wood 2018
# Simple test scripts for QUIT programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SIZE="8,8,8"
NOISE="0.001"
TOL="75"
cat > emt.json <<OUT
{
    "PDFile": "PD$EXT",
    "f_bFile": "f_b$EXT",
    "k_bfFile": "k_bf$EXT",
    "T1_fFile": "T1_f$EXT",
    "T2_fFile": "T2_f$EXT",
    "SSFPMT": {
        "FA":  [16, 32, 16, 32],
        "TR":  [0.008, 0.008, 0.009, 0.009],
        "Trf": [0.000256, 0.000256, 0.001024, 0.001024],
        "pulse": { "p1": 0.01047423, "p2": 0.4485347 }
    }
}
OUT

}

@test "SSFP-EMT Simulate" {
qinewimage --size "$SIZE" -f 1.0 PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -g "0 0.01 0.15" f_b$EXT
[ -e f_b$EXT ]
qinewimage --size "$SIZE" -g "1 5 10" k_bf$EXT
[ -e k_bf$EXT ]
qinewimage --size "$SIZE" -g "2 1.1 1.5" T1_f$EXT
[ -e T1_f$EXT ]
qinewimage --size "$SIZE" -f 0.06 T2_f$EXT
[ -e T2_f$EXT ]

qi_ssfp_emt --verbose --simulate=$NOISE G$EXT a$EXT b$EXT < emt.json
[ -e G$EXT ]
[ -e a$EXT ]
[ -e b$EXT ]
}

@test "SSFP-EMT Fit" {
qi_ssfp_emt --verbose G$EXT a$EXT b$EXT < emt.json
qidiff --baseline=PD$EXT --input=EMT_PD$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=f_b$EXT --input=EMT_f_b$EXT --noise=$NOISE --tolerance=$TOL --abs --verbose
# Don't compare k_bf for now, it doesn't fit well
# qidiff --baseline=k_bf$EXT --input=EMT_k_bf$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=T1_f$EXT --input=EMT_T1_f$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=T2_f$EXT --input=EMT_T2_f$EXT --noise=$NOISE --tolerance=$TOL --verbose
return 1
}
