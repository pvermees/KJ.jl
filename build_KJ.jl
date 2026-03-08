# 1. Define the app entry point
# This creates a temporary file that PackageCompiler will use to start your app.
entry_point_code = """
module KJApp
    using KJ
    function julia_main()::Cint
        try
            # By default, launch the TUI when the .exe is opened
            KJ.TUI()
        catch e
            Base.display_error(e)
            return 1 # Return non-zero on error
        end
        return 0 # Success
    end
end
"""

open("precompile_app.jl", "w") do f
    write(f, entry_point_code)
end

# 2. Set up paths based on the location of this script, not the current working directory
package_path = @__DIR__                           # Directory containing this build script
dest_path = joinpath(package_path, "KJ_Windows_Build")  # The output folder for your .exe

println("Building KJ.jl executable for Windows...")

# 3. Create a temporary precompile file and run the compilation
mktempdir() do tmpdir
    precompile_file = joinpath(tmpdir, "precompile_app.jl")
    open(precompile_file, "w") do f
        write(f, entry_point_code)
    end

    create_app(
        package_path,
        dest_path;
        # Points to the module containing julia_main
        executables = ["KJ" => "KJApp"],
        precompile_execution_file = precompile_file,
        force = true,
        filter_stdlibs = false,           # Set to true to reduce file size (~200MB saved)
        include_lazy_artifacts = true,    # Important if KJ uses specific data artifacts
        cpu_target = "generic;sandybridge,y;haswell,y" # Ensures compatibility across different Windows CPUs
    )
end

println("Build complete! Find your executable in: $dest_path/bin/KJ.exe")