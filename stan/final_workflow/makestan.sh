cd ~/cmdstan-2.15.0
module load GNU/6.2
make build
make O=3 ~/forestlight/stancode/model_ppow_withlik
make O=3 ~/forestlight/stancode/model_pexp_withlik
make O=3 ~/forestlight/stancode/model_wpow_withlik
make O=3 ~/forestlight/stancode/model_wexp_withlik
