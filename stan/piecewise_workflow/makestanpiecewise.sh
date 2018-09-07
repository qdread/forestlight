cd ~/cmdstan-2.15.0
module load GNU/6.2
make build
make O=3 ~/forestlight/stancode/model_d1p1
make O=3 ~/forestlight/stancode/model_d1p2
make O=3 ~/forestlight/stancode/model_d2p1
make O=3 ~/forestlight/stancode/model_d2p2
make O=3 ~/forestlight/stancode/model_d3p1
make O=3 ~/forestlight/stancode/model_d3p2
