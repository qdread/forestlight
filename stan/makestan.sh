cd ~/cmdstan-2.15.0
module load GNU/6.2
make build
make O=3 ~/forestlight/stancode/model_ppow
make O=3 ~/forestlight/stancode/model_pexp
make O=3 ~/forestlight/stancode/model_wpow
make O=3 ~/forestlight/stancode/model_wexp