# Copyright Tobias Wood 2017
# Tests for perfusion (ASL)

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
}

# @test "Perfusion (ASL)" {

# SIZE="32,32,32,2"
# qinewimage --verbose --dims=4 --size="$SIZE" --step="3 1 1.06 2" asl$EXT
# qi_asl --verbose asl$EXT <<END_INPUT
# {
#     "CASL" : {
#         "TR": 4.0,
#         "label_time": 3.0,
#         "post_label_delay": [0.3]
#     }
# }
# END_INPUT
# }

@test "Perfusion (ASE)" {
SIZE="3,3,3"
qinewimage S0$EXT --size="$SIZE" -f 1
qinewimage R2p$EXT --size="$SIZE" -g "0 1 5"
qinewimage DBV$EXT --size="$SIZE" -g "1 0.05 0.25"
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
    "DBVFile" : "DBV$EXT",
    "dTFile" : "dT$EXT"
}
END
SPIN_FILE="me$EXT"
NOISE="0.0"
qi_ase_oef --verbose --simulate=$NOISE $SPIN_FILE --threads=1 < input.json
return 1
qi_ase_oef --verbose $SPIN_FILE --threads=1 < input.json
qidiff --baseline=R2p$EXT --input=me_R2p$EXT --noise=$NOISE --tolerance=5 --verbose
qidiff --baseline=DBV$EXT --input=me_DBV$EXT --noise=$NOISE --tolerance=250 --abs --verbose
}

@test "Perfusion (Z-Shim)" {
# Pythagoras 3,4,5 triangle!
qinewimage zshim_in$EXT --dims=4 --size="2,2,2,2" --step="3 3 4 2"
qinewimage zshim_ref$EXT --size="2,2,2" --fill="5"
qi_zshim --verbose --zshims=2 zshim_in$EXT
qidiff --baseline=zshim_ref$EXT --input=zshim_in_zshim$EXT --abs --tolerance=0.1 --verbose
}