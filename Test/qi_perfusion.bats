# Copyright Tobias Wood 2017
# Tests for perfusion (ASL)

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
}

@test "Perfusion (ASL)" {

SIZE="32,32,32,2"
qinewimage --verbose --dims=4 --size="$SIZE" --step="3 1 1.06 2" asl$EXT
qi_asl --verbose asl$EXT <<END_INPUT
{
    "CASL" : {
        "TR": 4.0,
        "label_time": 3.0,
        "post_label_delay": [0.3]
    }
}
END_INPUT
}

@test "Perfusion (ASE)" {
SIZE="9,9,9"
qinewimage S0$EXT --size="$SIZE" -f 100
qinewimage R2p$EXT --size="$SIZE" -g "0 1 3"
qinewimage OEF$EXT --size="$SIZE" -g "1 0.25 0.5"
qinewimage dT$EXT --size="$SIZE" -g "2 -0.025 0.025"

cat > input.json <<END
{
    "MultiEchoFlex" : {
        "TR" : 2.5,
        "TE" : [ -0.05, 
                 -0.049,
                 -0.048,
                 -0.047,
                 -0.046,
                 -0.045,
                 -0.044,
                 -0.043,
                 -0.042,
                 -0.041,
                 -0.040,
                 -0.039,
                 -0.038,
                 -0.037,
                 -0.036,
                 -0.035,
                 -0.034,
                 -0.033,
                 -0.032,
                 -0.031,
                 -0.030,
                 -0.029,
                 -0.028,
                 -0.027,
                 -0.026,
                 -0.025,
                 -0.024,
                 -0.023,
                 -0.022,
                 -0.021,
                 -0.020,
                 -0.019,
                 -0.018,
                 -0.017,
                 -0.016,
                 -0.015,
                 -0.014,
                 -0.013,
                 -0.012,
                 -0.011,
                 -0.010,
                 -0.009,
                 -0.008,
                 -0.007,
                 -0.006,
                 -0.005,
                 -0.004,
                 -0.003,
                 -0.002,
                 -0.001,
                  0.000,
                  0.001,
                  0.002,
                  0.003,
                  0.004,
                  0.005,
                  0.006,
                  0.007,
                  0.008,
                  0.009,
                  0.010,
                  0.011,
                  0.012,
                  0.013,
                  0.014,
                  0.015,
                  0.016,
                  0.017,
                  0.018,
                  0.019,
                  0.020,
                  0.021,
                  0.022,
                  0.023,
                  0.024,
                  0.025,
                  0.026,
                  0.027,
                  0.028,
                  0.029,
                  0.030,
                  0.031,
                  0.032,
                  0.033,
                  0.034,
                  0.035,
                  0.036,
                  0.037,
                  0.038,
                  0.039,
                  0.040,
                  0.041,
                  0.042,
                  0.043,
                  0.044,
                  0.045,
                  0.046,
                  0.047,
                  0.048,
                  0.049,
                  0.050 ]
    },
    "S0File" : "S0$EXT",
    "R2pFile" : "R2p$EXT",
    "OEFFile" : "OEF$EXT",
    "dTFile" : "dT$EXT"
}
END
cat > input2.json <<END
{
    "MultiEchoFlex" : {
        "TR" : 2.5,
        "TE" : [ -0.05, 
                 -0.049,
                 -0.048,
                 -0.047,
                 -0.046,
                 -0.045,
                 -0.044,
                 -0.043,
                 -0.042,
                 -0.041,
                 -0.040,
                 -0.039,
                 -0.038,
                 -0.037,
                 -0.036,
                 -0.035,
                 -0.034,
                 -0.033,
                 -0.032,
                 -0.031,
                 -0.030,
                 -0.029,
                 -0.028,
                 -0.027,
                 -0.026,
                 -0.025,
                 -0.024,
                 -0.023,
                 -0.022,
                 -0.021,
                 -0.020,
                 -0.019,
                 -0.018,
                 -0.017,
                 -0.016,
                 -0.015,
                 -0.014,
                 -0.013,
                 -0.012,
                 -0.011,
                 -0.010,
                 -0.009,
                 -0.008,
                 -0.007,
                 -0.006,
                 -0.005,
                 -0.004,
                 -0.003,
                 -0.002,
                 -0.001,
                  0.000,
                  0.001,
                  0.002,
                  0.003,
                  0.004,
                  0.005,
                  0.006,
                  0.007,
                  0.008,
                  0.009,
                  0.010,
                  0.011,
                  0.012,
                  0.013,
                  0.014,
                  0.015,
                  0.016,
                  0.017,
                  0.018,
                  0.019,
                  0.020,
                  0.021,
                  0.022,
                  0.023,
                  0.024,
                  0.025,
                  0.026,
                  0.027,
                  0.028,
                  0.029,
                  0.030,
                  0.031,
                  0.032,
                  0.033,
                  0.034,
                  0.035,
                  0.036,
                  0.037,
                  0.038,
                  0.039,
                  0.040,
                  0.041,
                  0.042,
                  0.043,
                  0.044,
                  0.045,
                  0.046,
                  0.047,
                  0.048,
                  0.049,
                  0.050 ]
    },
    "S0File" : "me_S0$EXT",
    "R2pFile" : "me_R2p$EXT",
    "OEFFile" : "me_OEF$EXT",
    "dTFile" : "me_dT$EXT"
}
END
SPIN_FILE="me$EXT"
NOISE="0.1"
qi_ase_oef --verbose --simulate=$NOISE $SPIN_FILE --threads=1 < input.json
qi_ase_oef --verbose $SPIN_FILE --threads=1 < input.json
qi_ase_oef --verbose --simulate=$NOISE check.nii.gz --threads=1 < input2.json
qidiff --baseline=R2p$EXT --input=me_R2p$EXT --noise=$NOISE --tolerance=1 --verbose
qidiff --baseline=OEF$EXT --input=me_OEF$EXT --noise=$NOISE --tolerance=1 --verbose
}

@test "Perfusion (Z-Shim)" {
# Pythagoras 3,4,5 triangle!
qinewimage zshim_in$EXT --dims=4 --size="2,2,2,2" --step="3 3 4 2"
qinewimage zshim_ref$EXT --size="2,2,2" --fill="5"
qi_zshim --verbose --zshims=2 zshim_in$EXT
qidiff --baseline=zshim_ref$EXT --input=zshim_in_zshim$EXT --abs --tolerance=0.1 --verbose
}