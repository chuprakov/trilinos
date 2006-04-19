#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

my $base_base_dir = cwd();

my $base_options = 
" --moocho-options-file=$base_base_dir/Moocho.opt"
." --do-opt --use-first-order"
." --no-use-prec --lowsf=belos"
; # Can be overridden by --extra-args!

my $mpi_go = "mpirun -machinefile ~/bin/linux/machinelist-sahp6556 -np";

my $exe = "$base_base_dir/NLPThyraEpetraAdvDiffReactOpt.exe";

my $len_x = 1.0;
my $len_y = 1.0;
my $local_nx = 1;
my $global_ny = 2;

my $reaction_rate = 1.0;
my $beta = 1.0;

my $starting_num_procs = 1;
my $max_num_procs = 2;

my $starting_model_np = 1;
my $max_model_np = 1;
my $incr_model_np = 1;

my $starting_block_size = 1;
my $max_block_size = 1;
my $incr_block_size = 1;

my $extra_args = "";


GetOptions(
  "mpi-go=s"                => \$mpi_go,
  "exe=s"                   => \$exe,
  "len-x=i"                 => \$len_x,
  "len-y=i"                 => \$len_y,
  "local-nx=i"              => \$local_nx,
  "global-ny=i"             => \$global_ny,
  "reaction-rate=f"         => \$reaction_rate,
  "beta=f"                  => \$beta,
  "starting-num-procs=i"    => \$starting_num_procs,
  "max-num-procs=i"         => \$max_num_procs,
  "starting-model-np=i"     => \$starting_model_np,
  "max-model-np=i"          => \$max_model_np,
  "incr-model-np=i"         => \$incr_model_np,
  "starting-block-size=i"   => \$starting_block_size,
  "max-block-size=i"        => \$max_block_size,
  "incr-block-size=i"       => \$incr_block_size,
  "extra-args=s"            => \$extra_args
  );

my $base_dir_name = "runs";
mkdir $base_dir_name, 0777;
chdir $base_dir_name;

my $model_np = $starting_model_np;
my $vary_model_np = ( $model_np < $max_model_np );
for( ; $model_np <= $max_model_np; $model_np += $incr_model_np ) {

  if($vary_model_np) {
    my $dir_name = "model_np-$model_np";
    mkdir $dir_name, 0777;
    chdir $dir_name;
  }

  my $num_procs = $starting_num_procs;
  my $vary_num_procs = ( $num_procs < $max_num_procs );
  for( ; $num_procs <= $max_num_procs; $num_procs *= 2 ) {

    if($vary_num_procs) {
      my $dir_name = "num_procs-$num_procs";
      mkdir $dir_name, 0777;
      chdir $dir_name;
    }

    my $block_size = $starting_block_size;
    my $vary_block_size = ( $block_size < $max_block_size );
    for( ; $block_size <=  min($max_block_size,$model_np+1); $block_size += $incr_block_size ) {

      if($vary_block_size) {
        my $dir_name = "block_size-$block_size";
        mkdir $dir_name, 0777;
        chdir $dir_name;
      }

      my $curr_dir = cwd();

      print
        "\n***"
        ,"\n*** $curr_dir"
        ,"\n***\n";

      my $local_ny = $global_ny / $num_procs;

      my $cmnd =
        "$mpi_go $num_procs"
        ." $exe $base_options"
        ." --len-x=$len_x --len-y=$len_y --local-nx=$local_nx --local-ny=$local_ny"
        ." --np=$model_np"
        ." --reaction-rate=$reaction_rate --beta=$beta"
        ." --lowsf-extra-params=\"<ParameterList><Parameter name=\\\"Block Size\\\" type=\\\"int\\\" value=\\\"$block_size\\\"/></ParameterList>\""
        ." $extra_args"
        ." 2>&1 | tee run-test.out"
        ;

      run_cmnd($cmnd);

      if($vary_block_size) {
        chdir "..";
      }

    }

    if($vary_num_procs) {
      chdir "..";
    }

  }

  if($vary_model_np) {
    chdir "..";
  }

}

chdir "..";

#
# Subroutines
#

sub min {
  my $a = shift;
  my $b = shift;
  return ( $a < $b ? $a : $b );
}

sub run_cmnd {
    my  $cmnd = shift;
    print "\nRunning command: $cmnd\n";
    system($cmnd);
}
