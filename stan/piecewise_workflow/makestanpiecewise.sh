cd ~/cmdstan-2.15.0
# module load GNU/6.2 # no longer needed with new operating system
make build
make O=3 ~/forestlight/stancode/model_d1p1
make O=3 ~/forestlight/stancode/model_d1p2
make O=3 ~/forestlight/stancode/model_d2p1
make O=3 ~/forestlight/stancode/model_d2p2
make O=3 ~/forestlight/stancode/model_d3p1
make O=3 ~/forestlight/stancode/model_d3p2
make O=3 ~/forestlight/stancode/model_dwp1
make O=3 ~/forestlight/stancode/model_dwp2