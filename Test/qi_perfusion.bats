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
qinewimage dT$EXT --size="$SIZE" -g "0 -0.025 0.025"
qinewimage OEF$EXT --size="$SIZE" -g "1 0.2 0.5"
qinewimage DBV$EXT --size="$SIZE" -g "2 0.01 0.1"

cat > input.json <<END
{
    "MultiEchoFlex" : {
        "TR" : 2.5,
        "TE" : [ -0.05, 
                 -0.045,
                 -0.04,
                 -0.035,
                 -0.03,
                 -0.025,
                 -0.02,
                 -0.015,
                 -0.01,
                 -0.005,
                  0.0,
                  0.005,
                  0.01,
                  0.015,
                  0.02,
                  0.025,
                  0.03,
                  0.035,
                  0.04,
                  0.045,
                  0.05 ]
    },
    "S0File" : "S0$EXT",
    "dTFile" : "dT$EXT",
    "OEFFile" : "OEF$EXT",
    "DBVFile" : "DBV$EXT"
}
END
cat > input2.json <<END
{
    "MultiEchoFlex" : {
        "TR" : 2.5,
        "TE" : [ -0.05, 
                 -0.045,
                 -0.04,
                 -0.035,
                 -0.03,
                 -0.025,
                 -0.02,
                 -0.015,
                 -0.01,
                 -0.005,
                  0.0,
                  0.005,
                  0.01,
                  0.015,
                  0.02,
                  0.025,
                  0.03,
                  0.035,
                  0.04,
                  0.045,
                  0.05 ]
    },
    "S0File" : "me_S0$EXT",
    "dTFile" : "me_dT$EXT",
    "OEFFile" : "me_OEF$EXT",
    "DBVFile" : "me_DBV$EXT"
}
END
SPIN_FILE="me$EXT"
NOISE="0.01"
qi_ase_oef --verbose --simulate=$NOISE $SPIN_FILE --threads=1 < input.json
qi_ase_oef --verbose $SPIN_FILE --threads=1 < input.json
qi_ase_oef --verbose --simulate=$NOISE check.nii.gz --threads=1 < input2.json
qidiff --baseline=OEF$EXT --input=me_OEF$EXT --noise=$NOISE --tolerance=5 --verbose
qidiff --baseline=DBV$EXT --input=me_DBV$EXT --noise=$NOISE --tolerance=5 --verbose
}

@test "Perfusion (Z-Shim)" {
# Pythagoras 3,4,5 triangle!
qinewimage zshim_in$EXT --dims=4 --size="2,2,2,2" --step="3 3 4 2"
qinewimage zshim_ref$EXT --size="2,2,2" --fill="5"
qi_zshim --verbose --zshims=2 zshim_in$EXT
qidiff --baseline=zshim_ref$EXT --input=zshim_in_zshim$EXT --abs --tolerance=0.1 --verbose
}