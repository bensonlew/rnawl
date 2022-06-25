import os
def envir(bind_obj, type="perl"):
    base_path = bind_obj.config.SOFTWARE_DIR
    if type not in ["perl", "python", "magic", "c", "r"]:
        bind_obj.logger.info("TAKE CARE, no environmental config for type: %s" % type)
    if type == "perl":
        bind_obj.set_environ(PERL5LIB="/mnt/lustre/users/sanger/app/program/perl-5.24.0/lib/site_perl/5.24.0:/mnt/lustre/users/sanger/app/program/perl/perls/perl-5.24.0/lib:/mnt/lustre/users/sanger/app/program/perl/perls/perl-5.24.0/lib/5.24.0/x86_64-linux-thread-multi:/mnt/lustre/users/sanger/app/program/perl/perls/perl-5.24.0/lib/site_perl/5.24.0/x86_64-linux-thread-multi:/mnt/lustre/users/sanger/app/program/perl/perls/perl-5.24.0/lib")
    elif type == "python":
        pass
    elif type == "magic":
        pass
    elif type == "c":
        pass
    elif type == "r":
        # bind_obj.set_environ(LD_LIBRARY_PATH="/mnt/lustre/users/sanger/app/library/lib")
        bind_obj.set_environ(PATH="/mnt/lustre/users/sanger/app/program/R-3.3.1/bin", LD_LIBRARY_PATH="/mnt/lustre/users/sanger/app/program/R-3.3.1/lib64/R/lib", LIBRARY_PATH="/mnt/lustre/users/sanger/app/program/R-3.3.1/lib64/R/lib", C_INCLUDE_PATH="/mnt/lustre/users/sanger/app/program/R-3.3.1/lib64/R/include", CPLUS_INCLUDE_PATH="/mnt/lustre/users/sanger/app/program/R-3.3.1/lib64/R/include")
    elif type == "bash_profile":
        path_list = [
            os.path.join(base_path, "program/nano"),
            os.path.join(base_path, "program/sun_jdk1.8.0/bin"),
            os.path.join(base_path, "program/R-3.3.1/bin"),
            os.path.join(base_path, "miniconda2/bin"),
            os.path.join(base_path, "program/vim/bin"),
            "/mnt/lustre/users/sanger/biocluster/bin",
            os.path.join(base_path, "library/bin")
        ]
        ld_library_path_list = [
            os.path.join(base_path, "program/R-3.3.1/lib64/R/lib"),
            os.path.join(base_path, "bioinfo/seq/EMBOSS-6.6.0/lib:/opt/openmpi/lib")
        ]
        library_path_list = [
            os.path.join(base_path, "program/R-3.3.1/lib64/R/lib"),
            os.path.join(base_path, "library/lib"),
            os.path.join(base_path, "library/lib64")
        ]
        c_include_path_list = [
            os.path.join(base_path, "program/R-3.3.1/lib64/R/include"),
            os.path.join(base_path, "program/Python/include"),
            os.path.join(base_path, "library/include")
        ]
        bind_obj.set_environ(PATH=":".join(path_list),
                             LD_LIBRARY_PATH=":".join(ld_library_path_list),
                             LIBRARY_PATH=":".join(library_path_list),
                             CPLUS_INCLUDE_PATH=":".join(c_include_path_list),
                             C_INCLUDE_PATH=":".join(c_include_path_list))