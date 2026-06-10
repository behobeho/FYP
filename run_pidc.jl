input_dir = raw"C:\Users\bella\OneDrive - Imperial College London\FYP\Code\BENCHMARKING\HVG_filtered"
output_dir = raw"C:\Users\bella\OneDrive - Imperial College London\FYP\Code\BENCHMARKING\Final_Filtered"

for i in 0:13

	println("Processing file $i")

	input_file = joinpath(input_dir, "100_HVG_$(i).csv")

	df = CSV.read(input_file, DataFrame)
	df = DataFrame(permutedims(Matrix(df)), :auto)
		
	temp_file = "temp_PIDC_input.csv"
	CSV.write(temp_file, df)

	nodes = get_nodes(temp_file, delim=',')

	@time inferred_network = InferredNetwork(
		PIDCNetworkInference(),
		nodes
	)

	A = get_adjacency_matrix(inferred_network, 0.1)

	adj_matrix, name_to_idx, idx_to_name = A

	n = size(adj_matrix, 1)

	gene_names = [idx_to_name[j] for j in 1:n]

	result_df = DataFrame(adj_matrix .* 1.0, Symbol.(gene_names))

	insertcols!(result_df, 1, :gene => gene_names)

	# Output CSV
	output_file = joinpath(
		output_dir,
		"100_HVG_$(i)_PIDC.csv"
	)

	CSV.write(output_file, result_df)

	println("Saved to $output_file")

end