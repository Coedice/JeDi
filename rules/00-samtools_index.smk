#######################################################################################
rule samtools_index:
	input:
		config['gstacks']['input_dir'] + '{xyz}.bam'
	output:
		config['gstacks']['input_dir'] + '{xyz}.bam.csi'
	threads:
		config['threads']
	shell:
		"samtools index -c -@ {threads} {input}"
