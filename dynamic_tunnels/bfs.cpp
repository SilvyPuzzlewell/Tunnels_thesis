//Modified version of http://www.geeksforgeeks.org/breadth-first-traversal-for-a-graph/






bool search(shared_ptr<vertex> base, shared_ptr<vertex> target){

    int frame_diference = find_frame_difference;
    

    double  direction_vector_length = compute_metric_eucleidean(target->get_location_coordinates(), base->get_location_coordinates());
    int num_of_interpolated_steps = direction_vector_length / min_step;

    if(direction_vector_length <= min_step){
        return false;
    }

    double* increment_vector = create_direction_vector(base, target, min_step);
    base->get_direction_coordinates();

    vector<vector<unsigned char>> matrix(frame_diference, vector<int>(num_of_interpolated_steps));

    for(int i = 0; i < num_of_interpolated_steps; i++){
        double* current_coordinates = add_vector(base->get_location_coordinates, increment_vector, ADDITION);
        for(int j = 0; j <= frame_diference; j++){
            if(is_in_obstacle_custom_frame(current_coordinates, probe_radius, DONT_CHECK_WITH_BLOCKING_SPHERES, base->get_first_frame() + j)){
                matrix[i][j] = 0;
            }
            else {
                matrix[i][j] = 1;
            } 
        }
        delete [] current_coordinates;
    }

    queue <array<int, 2>> visited;
    array<int, 2> cur = {0, 0};
    visited.push(queue);

    bool success = false;

    while(visited.size() != 0){
        if(cur[0] < num_of_interpolated_steps - 1){
            if(matrix[cur[0] + 1][cur[1]] == 1){
                if(cur[0] == num_of_interpolated_steps - 1 && cur[1] = frame_diference - 1){
                    success = true;
                    break;
                }
                array<int, 2> new_position = {cur[0] + 1, cur[1]}
                visited.push(new_position);
            }
        }   

        if(cur[0] < frame_diference - 1){
            if(matrix[cur[0]][cur[1] + 1] == 1){
                if(cur[0] == num_of_interpolated_steps - 1 && cur[1] = frame_diference - 1){
                    success = true;
                    break;
                }
                array<int, 2> new_position = {cur[0], cur[1] + 1}
                visited.push(new_position);
            }
        }

        cur = visited.front();
        visited.pop();

    }
    return success;
}




