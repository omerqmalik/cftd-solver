join -j1 -o1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.4,2.5,2.6,2.8,2.7 -t ','  \
	<(<outside_log_file.txt awk -F "," '{print $1"-"$2","$0}' | sort -k1,1) \
	<(<inside_mem_log.txt awk -F "," '{print $1"-"$2","$0}' | sort -k1,1)