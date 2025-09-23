"""
CAT-Surface Python Wrapper Examples

This file contains comprehensive examples of how to use the CAT-Surface Python wrapper
for various neuroimaging surface analysis tasks.
"""

import os
from cat_surface import (
    set_executable_path, CATConfig, CATSurfaceError,
    # Volume processing
    segment_with_amap, denoise_sanlm, compute_thickness_pbt,
    extract_surface_marching_cubes, project_volume_to_surface,
    # Surface processing
    average_surfaces, create_pial_white_surfaces, calculate_surface_distance,
    warp_surface, create_spherical_surface, resample_surface
)
from cat_surface.volume import VolAmap, VolThicknessPbt
from cat_surface.surface import SurfAverage, SurfDistance


def setup_cat_surface():
    """
    Setup CAT-Surface Python wrapper with custom executable path.
    """
    # Option 1: Set global executable path
    set_executable_path("/Users/gaser/Library/CloudStorage/Dropbox/GitHub/CAT-Surface/build-native-arm64/Progs")
    
    # Option 2: Create custom configuration
    config = CATConfig("/Users/gaser/Library/CloudStorage/Dropbox/GitHub/CAT-Surface/build-native-arm64/Progs")
    return config


def example_volume_segmentation():
    """
    Example: Segment brain volume using adaptive MAP (CAT_VolAmap).
    """
    print("=== Volume Segmentation Example ===")
    
    # Input files (replace with your actual file paths)
    input_volume = "data/brain_t1.nii"
    label_file = "data/initial_labels.nii"
    output_basename = "output/segmented"
    
    try:
        # Method 1: Using convenience function
        result = segment_with_amap(
            input_file=input_volume,
            label_file=label_file,
            output_file=output_basename,
            iters=100,                    # Increase iterations for better accuracy
            pve=True,                     # Use partial volume estimation
            cleanup=2,                    # Strong cleanup
            write_seg=(True, True, True), # Write CSF, GM, WM separately
            verbose=True
        )
        
        print(f"Segmentation completed successfully!")
        print(f"Return code: {result.returncode}")
        print(f"Output files: {result.output_files}")
        
        # Method 2: Using class interface for more control
        segmenter = VolAmap()
        result2 = segmenter.segment(
            input_file=input_volume,
            label_file=label_file,
            output_file=f"{output_basename}_detailed",
            mrf_weight=0.1,       # Add MRF regularization
            bias_fwhm=30.0,       # Apply bias correction
            write_bias=True       # Save bias field
        )
        
    except CATSurfaceError as e:
        print(f"Error during segmentation: {e}")
    except FileNotFoundError as e:
        print(f"Input file not found: {e}")


def example_thickness_estimation():
    """
    Example: Compute cortical thickness using projection-based method.
    """
    print("=== Cortical Thickness Estimation Example ===")
    
    # Input files
    pve_labels = "data/pve_labels.nii"
    gmt_output = "output/thickness.nii"
    ppm_output = "output/position_map.nii"
    
    try:
        # Compute thickness with custom parameters
        result = compute_thickness_pbt(
            input_file=pve_labels,
            gmt_output=gmt_output,
            ppm_output=ppm_output,
            n_avgs=10,                # More averaging for smoother results
            fill_holes=0.3,           # Threshold for hole filling
            correct_voxelsize=0.5,    # Correct for voxel size effects
            verbose=True
        )
        
        print(f"Thickness estimation completed!")
        print(f"Thickness map: {gmt_output}")
        print(f"Position map: {ppm_output}")
        
        # Also compute additional outputs
        thickness_tool = VolThicknessPbt()
        result_extended = thickness_tool.compute_thickness(
            input_file=pve_labels,
            gmt_output="output/thickness_ext.nii",
            ppm_output="output/position_ext.nii",
            wmd_output="output/wm_distance.nii",    # White matter distance
            csd_output="output/csf_distance.nii",   # CSF distance
            sharpen=0.1                             # Sharpen PPM map
        )
        
    except CATSurfaceError as e:
        print(f"Error during thickness computation: {e}")


def example_surface_extraction():
    """
    Example: Extract cortical surface using marching cubes.
    """
    print("=== Surface Extraction Example ===")
    
    # Extract GM/WM boundary surface
    pve_volume = "data/pve_labels.nii"
    surface_output = "output/central_surface.gii"
    
    try:
        result = extract_surface_marching_cubes(
            input_file=pve_volume,
            output_file=surface_output,
            threshold=1.5  # GM/WM boundary threshold
        )
        
        print(f"Surface extracted: {surface_output}")
        
    except CATSurfaceError as e:
        print(f"Error during surface extraction: {e}")


def example_surface_processing():
    """
    Example: Process cortical surfaces - create pial/white, calculate distances.
    """
    print("=== Surface Processing Example ===")
    
    # Input files
    central_surface = "data/central.gii"
    thickness_values = "data/thickness.txt"
    
    try:
        # Create pial and white matter surfaces
        result = create_pial_white_surfaces(
            central_surface=central_surface,
            thickness_file=thickness_values,
            pial_output="output/pial.gii",
            white_output="output/white.gii",
            weight=0.5,           # Use half thickness for displacement
            equivolume=True,      # Apply equi-volume correction
            verbose=True
        )
        
        print("Pial and white surfaces created!")
        
        # Calculate distance between surfaces
        distance_result = calculate_surface_distance(
            surface1="output/white.gii",
            surface2="output/pial.gii",
            output_values="output/surface_distance.txt"
        )
        
        print("Surface distance calculated!")
        
        # Alternative: Use thickness file directly
        distance_from_thickness = calculate_surface_distance(
            surface1=central_surface,
            surface2="",  # Not used when thickness_file is provided
            output_values="output/thickness_distance.txt",
            thickness_file=thickness_values
        )
        
    except CATSurfaceError as e:
        print(f"Error during surface processing: {e}")


def example_surface_analysis():
    """
    Example: Advanced surface analysis - averaging, spherical mapping, resampling.
    """
    print("=== Advanced Surface Analysis Example ===")
    
    try:
        # Average multiple surfaces
        surface_list = [
            "data/subject1_surface.gii",
            "data/subject2_surface.gii",
            "data/subject3_surface.gii"
        ]
        
        avg_result = average_surfaces(
            surface_files=surface_list,
            output_file="output/average_surface.gii"
        )
        
        print("Surface averaging completed!")
        
        # Create spherical parameterization
        sphere_result = create_spherical_surface(
            input_surface="output/average_surface.gii",
            output_sphere="output/sphere.gii",
            verbose=True
        )
        
        print("Spherical parameterization created!")
        
        # Resample surface to template
        resample_result = resample_surface(
            input_surface="data/individual_surface.gii",
            template_surface="templates/standard_surface.gii",
            output_surface="output/resampled_surface.gii",
            values_file="data/surface_values.txt",
            output_values="output/resampled_values.txt"
        )
        
        print("Surface resampling completed!")
        
    except CATSurfaceError as e:
        print(f"Error during surface analysis: {e}")


def example_batch_processing():
    """
    Example: Batch processing multiple subjects.
    """
    print("=== Batch Processing Example ===")
    
    subjects = ["subj001", "subj002", "subj003"]
    
    for subject in subjects:
        print(f"Processing {subject}...")
        
        try:
            # Define subject-specific paths
            input_vol = f"data/{subject}/t1.nii"
            labels = f"data/{subject}/labels.nii"
            output_dir = f"output/{subject}"
            
            # Ensure output directory exists
            os.makedirs(output_dir, exist_ok=True)
            
            # Step 1: Segmentation
            seg_result = segment_with_amap(
                input_file=input_vol,
                label_file=labels,
                output_file=f"{output_dir}/segmented",
                verbose=False  # Reduce output for batch processing
            )
            
            # Step 2: Thickness estimation
            thickness_result = compute_thickness_pbt(
                input_file=f"{output_dir}/segmented_label.nii",
                gmt_output=f"{output_dir}/thickness.nii",
                ppm_output=f"{output_dir}/position.nii"
            )
            
            # Step 3: Surface extraction
            surface_result = extract_surface_marching_cubes(
                input_file=f"{output_dir}/segmented_label.nii",
                output_file=f"{output_dir}/surface.gii"
            )
            
            print(f"  {subject} completed successfully!")
            
        except CATSurfaceError as e:
            print(f"  Error processing {subject}: {e}")
            continue


def example_error_handling():
    """
    Example: Proper error handling and debugging.
    """
    print("=== Error Handling Example ===")
    
    try:
        # This will fail - file doesn't exist
        result = segment_with_amap(
            input_file="nonexistent_file.nii",
            label_file="also_nonexistent.nii"
        )
        
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        
    except CATSurfaceError as e:
        print(f"CAT-Surface error: {e}")
        if hasattr(e, 'stderr') and e.stderr:
            print(f"Error details: {e.stderr}")
            
    # Using class interface for more detailed error handling
    segmenter = VolAmap()
    try:
        result = segmenter.segment(
            input_file="data/brain.nii",
            label_file="data/labels.nii",
            iters=-1  # Invalid parameter
        )
    except CATSurfaceError as e:
        print(f"Process error: {e}")
        print(f"Return code: {e.returncode}")
        print(f"Command: {' '.join(e.cmd)}")


def example_volume_to_surface_projection():
    """
    Example: Project volume data onto surface.
    """
    print("=== Volume to Surface Projection Example ===")
    
    try:
        # Project functional data onto anatomical surface
        result = project_volume_to_surface(
            volume_file="data/functional_activation.nii",
            surface_file="data/central_surface.gii",
            output_file="output/surface_activation.txt"
        )
        
        print("Volume to surface projection completed!")
        
    except CATSurfaceError as e:
        print(f"Error during projection: {e}")


def main():
    """
    Run all examples (comment out sections you don't want to run).
    """
    print("CAT-Surface Python Wrapper Examples")
    print("=" * 40)
    
    # Setup
    try:
        config = setup_cat_surface()
        print("CAT-Surface configuration successful!")
    except Exception as e:
        print(f"Setup failed: {e}")
        print("Please adjust the executable path in setup_cat_surface()")
        return
    
    # Run examples (comment out as needed)
    # example_volume_segmentation()
    # example_thickness_estimation()
    # example_surface_extraction()
    # example_surface_processing()
    # example_surface_analysis()
    # example_batch_processing()
    example_error_handling()  # Safe to run without actual data
    # example_volume_to_surface_projection()
    
    print("\\nExamples completed!")


if __name__ == "__main__":
    main()