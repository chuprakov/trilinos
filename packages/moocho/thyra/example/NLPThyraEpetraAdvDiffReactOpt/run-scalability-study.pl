#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

my $len_x                  = 1.0;
my $local_len_y            = 1.0;
my $local_nx               = 4;
my $local_ny               = 4;
my $reaction_rate          = 1e-3;
my $beta                   = 0.0;
my $max_linear_iters       = 10;
my $max_rsqp_iters         = 1;
my $starting_num_procs     = 1;
my $max_num_procs          = 1;
my $exe                    = "";
my $mpi_go                 = "";
my $has_expat              = 1;
my $extra_args             = "";

GetOptions(
  "len-x=i"                 => \$len_x,
  "local-len-y=i"           => \$local_len_y,
  "local-nx=i"              => \$local_nx,
  "local-ny=i"              => \$local_ny,
  "reaction-rate=f"         => \$reaction_rate,
  "beta=f"                  => \$beta,
  "max-linear-iters=i"      => \$max_linear_iters,
  "max-rsqp-iters=i"        => \$max_rsqp_iters,
  "starting-num-procs=i"    => \$starting_num_procs,
  "max-num-procs=i"         => \$max_num_procs,
  "exe=s"                   => \$exe,
  "mpi-go=s"                => \$mpi_go,
  "has-expat!"              => \$has_expat,
  "extra-args=s"            => \$extra_args
  );

my $base_base_dir = cwd();

if($exe eq "") {
    $exe = "$base_base_dir/NLPThyraEpetraAdvDiffReactOpt.exe";
}

if($mpi_go eq "") {
    $mpi_go = "mpirun -np";
}

my $base_dir_name = "runs";
mkdir $base_dir_name, 0777;
chdir $base_dir_name;

my $num_procs = $starting_num_procs;

for( ; $num_procs <= $max_num_procs; $num_procs *= 2 ) {

    my $dir_name = "num_procs-$num_procs";
    mkdir $dir_name, 0777;
    chdir $dir_name;

    my $curr_dir = cwd();

    print
        "\n***"
        ,"\n*** $curr_dir"
        ,"\n***\n";

    my $len_y = $local_len_y * $num_procs;
    my $belosParams = getBelosParams();
    my $moochoOptions = getMoochoOptions();

    my $cmnd =
        "$mpi_go $num_procs"
        ." $exe"
        ." --no-use-prec --do-opt --use-direct --lowsf=belos"
        .( $has_expat ? " --lowsf-extra-params=\"$belosParams\"" : "" )
        ." --len-x=$len_x --len-y=$len_y --local-nx=$local_nx --local-ny=$local_ny"
        ." --reaction-rate=$reaction_rate --beta=$beta"
        ." --x0=0.0 --np=1 --p0=0.0"
        ." --moocho-extra-options=\"$moochoOptions\""
        ." $extra_args"
        ." | tee run-test.out"
        ;
    
    run_cmnd($cmnd);

    chdir "..";

}

chdir "..";

#
# Subroutines
#

sub getBelosParams {
    return
        "<ParameterList>"
        ."  <Parameter name=\\\"Default Rel Res Norm\\\" type=\\\"double\\\" value=\\\"1e-20\\\"/>"
        ."  <!-- Above: Set small tolerance to make sure the hit max iters. -->"
        ."  <Parameter name=\\\"Max Iters\\\" type=\\\"int\\\" value=\\\"$max_linear_iters\\\"/>"
        ."  <Parameter name=\\\"Block Size\\\" type=\\\"int\\\" value=\\\"1\\\"/>"
        ."  <ParameterList name=\\\"GMRES\\\">"
        ."    <Parameter name=\\\"Max Number of Krylov Vectors\\\" type=\\\"int\\\" value=\\\"$max_linear_iters\\\"/>"
        ."  </ParameterList>"
        ."  <ParameterList name=\\\"Outputter\\\">"
        ."    <Parameter name=\\\"Output Frequency\\\" type=\\\"int\\\" value=\\\"10\\\"/>"
        ."    <Parameter name=\\\"Output Max Res Only\\\" type=\\\"bool\\\" value=\\\"1\\\"/>"
        ."  </ParameterList>"
        ."</ParameterList>"
        ;
}

sub getMoochoOptions {
    return
        "MoochoSolver{test_nlp=false}"
        .":NLPSolverClientInterface{"
        ."opt_tol=1e-20,max_iter=$max_rsqp_iters,max_run_time=20"
        .",journal_output_level=PRINT_ALGORITHM_STEPS"
        .",calc_conditioning=true"
        .",calc_matrix_info_null_space_only=true"
        ."}"
        .":IterationPack::Algorithm{interrupt_file_name=\\\"interrupt.in\\\"}"
        .":DecompositionSystemStateStepBuilderStd{null_space_matrix=EXPLICIT,range_space_matrix=ORTHOGONAL}"
        .":NLPAlgoConfigMamaJama{line_search_method=FILTER}"
        .":LineSearchFilter{f_min=0.0}"
        ;
}

sub run_cmnd {
    my  $cmnd = shift;
    print "\nRunning command: $cmnd\n";
    system($cmnd);
}
