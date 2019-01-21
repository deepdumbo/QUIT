# Copyright Tobias Wood 2018
# Tests for various utilities

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "RF Profile" {

SIZE="32, 32, 32"
qinewimage --size "$SIZE" --grad="1 0 1" grad$EXT
cat > rf.json <<END
{
    "rf_pos" :  [0, 1],
    "rf_vals" : [0, 1]
}
END
qi_rfprofile grad$EXT b1$EXT --verbose --file=rf.json
}

@test "K-Space Filtering" {

SIZE="64,64,16,4"
qinewimage --dims=4 --size="$SIZE" --step="0 0 8 4" steps$EXT
qikfilter steps$EXT --threads=1 --filter_per_volume --filter=Gauss,2.0 --filter=Blackman --filter=Hamming --filter=Tukey --verbose
}