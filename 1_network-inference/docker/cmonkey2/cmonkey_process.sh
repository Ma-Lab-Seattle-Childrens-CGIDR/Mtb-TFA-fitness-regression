#!/bin/bash

db_file="$1"
tf_file="$2"
out_file="$3"
run="$(date +%s)"

out_dir="$(dirname $out_file)"

echo "Extracting data from sqlite database $db_file..."

read -d '\n' sql_query << EOT
select
	row_names.name,
	row_members.cluster,
	row_members.iteration,
	cluster_stats.residual
from
	row_members
inner join
	row_names
		using(order_num)
inner join
	cluster_stats
		on row_members.cluster = cluster_stats.cluster
		and row_members.iteration = cluster_stats.iteration;
EOT

sqlite3 -csv "$db_file" "$sql_query" > "$out_dir/cmonkey2-raw-${run}.csv"

echo "Done."
echo "Aggregating results across iterations..."

python3 /cmonkey/cmonkey_process.py \
	"$out_dir/cmonkey2-raw-${run}.csv" \
	"$tf_file" \
	"$out_file"

rm "$out_dir/cmonkey2-raw-${run}.csv"

echo "Done. Results saved to $out_file."
