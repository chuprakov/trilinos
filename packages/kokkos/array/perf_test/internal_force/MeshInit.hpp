template<class scalar_mdarray, class int_mdarray>
void mesh_init(	scalar_mdarray 	& coords, 
				int_mdarray 	& nodeIDs, 
				int_mdarray 	& elemIDs, 
				int_mdarray 	& epn_cumulative_counter,
				int elems_x, int elems_y, int elems_z){

//	construct a structured, rectangular prism mesh of Hex elements, 
//	with dimensions given by elems_x, elems_y, elems_z

	int nx = elems_x + 1;
	int ny = elems_y + 1;
	int nz = elems_z + 1;
	int nnodes = nx * ny * nz;
	int nelems = elems_x * elems_y * elems_z;

	nodeIDs = Kokkos::create_mdarray< int_mdarray >(nelems, 8);
	coords  = Kokkos::create_mdarray< scalar_mdarray >(nelems, 3, 8);
	epn_cumulative_counter = Kokkos::create_mdarray< int_mdarray >(nnodes + 1);

	int_mdarray elem_per_node_counter = Kokkos::create_mdarray< int_mdarray >(nnodes);

	for(int i = 0; i < elems_x; i++){
		for(int j = 0; j < elems_y; j++){
			for(int k = 0; k < elems_z; k++){

				int adj_index  = k * nx * ny + j * nx + i;
				int elem_index = k * elems_x * elems_y + j * elems_x + i;

				coords(elem_index, 0, 0) = i;
				coords(elem_index, 1, 0) = j;				
				coords(elem_index, 2, 0) = k;

				coords(elem_index, 0, 1) = i + 1;
				coords(elem_index, 1, 1) = j;				
				coords(elem_index, 2, 1) = k;

				coords(elem_index, 0, 2) = i + 1;
				coords(elem_index, 1, 2) = j;				
				coords(elem_index, 2, 2) = k - 1;

				coords(elem_index, 0, 3) = i;
				coords(elem_index, 1, 3) = j;				
				coords(elem_index, 2, 3) = k - 1;

				coords(elem_index, 0, 4) = i;
				coords(elem_index, 1, 4) = j + 1;				
				coords(elem_index, 2, 4) = k;

				coords(elem_index, 0, 5) = i + 1;
				coords(elem_index, 1, 5) = j + 1;				
				coords(elem_index, 2, 5) = k;

				coords(elem_index, 0, 6) = i + 1;
				coords(elem_index, 1, 6) = j + 1;				
				coords(elem_index, 2, 6) = k - 1;

				coords(elem_index, 0, 7) = i;
				coords(elem_index, 1, 7) = j + 1;				
				coords(elem_index, 2, 7) = k - 1;

				
				nodeIDs(elem_index, 0) = adj_index;
				nodeIDs(elem_index, 1) = adj_index 						+ 1;
				nodeIDs(elem_index, 2) = adj_index 			+ nx * ny 	+ 1;
				nodeIDs(elem_index, 3) = adj_index 			+ nx * ny;
				nodeIDs(elem_index, 4) = adj_index + nx;
				nodeIDs(elem_index, 5) = adj_index + nx 				+ 1;
				nodeIDs(elem_index, 6) = adj_index + nx 	+ nx * ny 	+ 1;
				nodeIDs(elem_index, 7) = adj_index + nx 	+ nx * ny;

			//	in addition to storing which nodes belong to each element,
			//	construct the transposed data structure: which elements a
			//	given node is in.

			//	To later gather directly into a CSR matrix, we start by 
			//	counting how many elements each node belongs to

				elem_per_node_counter(adj_index)++;
				elem_per_node_counter(adj_index + 1)++;
				elem_per_node_counter(adj_index + nx)++;
				elem_per_node_counter(adj_index + nx + 1)++;
				elem_per_node_counter(adj_index + nx * ny)++;
				elem_per_node_counter(adj_index + nx * ny + 1)++;
				elem_per_node_counter(adj_index + nx * ny + nx)++;
				elem_per_node_counter(adj_index + nx * ny + nx + 1)++;

			}
		}
	}

	int sum = 0;
	epn_cumulative_counter(0) = 0;

	for(int i = 0; i < nnodes; i++){

	//	Then, translate the counter to a cumulative one
	//	to set up the elemIDs data structure
		sum += elem_per_node_counter(i);		
		epn_cumulative_counter(i + 1) = sum;
		elem_per_node_counter(i) = 0;
		
	}

	elemIDs = Kokkos::create_mdarray< int_mdarray >(sum, 2);

//	for each element..
	for(int i = 0; i < elems_x; i++){
		for(int j = 0; j < elems_y; j++){
			for(int k = 0; k < elems_z; k++){

				int elem_index = k * elems_x * elems_y + j * elems_x + i;

			//	scatter its element ID information into the elemIDs structure
				for(int n = 0; n < 8; n++){

					int nid = nodeIDs(elem_index, n);

					elemIDs(epn_cumulative_counter(nid) + elem_per_node_counter(nid), 0) = elem_index;
					elemIDs(epn_cumulative_counter(nid) + elem_per_node_counter(nid), 1) = n;
					elem_per_node_counter(nid)++;
					
				}				

			}
		}
	}

}
