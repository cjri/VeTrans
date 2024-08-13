awk 'NR>2 {print last} {last=$0}' $1 | cut -f1 -d " " > Previous_haps.in
