#!/usr/bin/env perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use Storable;
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fKey,$genotype,$dOut,$start,$end,$stepSize,$minGroup,$maxGroup,$nChro,$redo,$log,$mode,$dumperboo);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"g:s"=>\$genotype,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				"n:s"=>\$nChro,

				"b:s"=>\$start,
				"e:s"=>\$end,
				"s:s"=>\$stepSize,
				"minGroup:s"=>\$minGroup,
				"maxGroup:s"=>\$maxGroup,

				"redo:i"=>\$redo,
				"mode:i"=>\$mode,
				"dumper:i"=>\$dumperboo,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $nChro);
#-------------------------------------------------------------------
#Global parameter
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

$start||=3;        ## lod start 
$end||=20;         ## lod end 
$stepSize||=1;     ## lod step size
$minGroup||=20;    ## minimum threshold for group size
$maxGroup||=500;   ## maximum threshold for group size
$|=1;              ## close buffer 
$redo //=0;##<
$mode||=2;         ## choose grouping scheme with minimum variance or maximum marker numbers: 1: max number; 2��min variance 
$dumperboo//=0;     ## choose whether print tree hashes to dumper 

if ($minGroup > $maxGroup) {	
	print "Illegal parameters' value for minGroup and maxGroup,minGroup > maxGroup.\n";
	exit;
}

open ($log,">$dOut/$fKey.linkageGroupping.log") or die $!;
print_log_head($log);
print $log "++Now lod=($start,$end),redo=$redo++++\n";

#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my %mLOD = ();       ## hash to store pair wise data
#my %Tree = ();       ## Group tree 
my $Tree = {};        ## group tree hash 
my %all_loci = ();   ## all loci set 
my $succeed_flag = 0; ## is last step finished 

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------

#my %info;
if (defined $genotype) {
    open IN,$genotype;
	while (<IN>) {
		chomp;
		next if (/^$/ || /MarkerID/ || /^name/ || /^popt/ || /^nloc/ || /^nind/ || /^#/);
		my ($marker) = split;
		$all_loci{$marker} = 1;
	}
    close IN;
}
#print Dumper %all_loci;die;

open (IN,"$fIn") or die $!;
<IN>;
while (<IN>) {
	chomp;
	next if (/^$/ or /^;/ || /^#/) ;
	my ($loci_1,$loci_2,$mLOD) = split;
	if (defined $genotype) {
		next if (!exists $all_loci{$loci_1} || !exists $all_loci{$loci_2}); ## filter the info of marker not found in genotype file 
	}
	$mLOD{$loci_1}{$loci_2} = $mLOD;
	$mLOD{$loci_2}{$loci_1} = $mLOD;
#	print "++$mLOD{$loci_1}{$loci_2}++\n";
	if (not defined $genotype){
		$all_loci{$loci_1} = 1;
		$all_loci{$loci_2} = 1;
	}
}
close IN;
#$open DUMP,">$dOut/$fKey.mloddumped";
#$print DUMP Dumper %mLOD;
#$close DUMP;

#-------------------------------------------------------------------
# determine linkage group 
#-------------------------------------------------------------------

#####################################################################
#                                                                   #
#                CONSTRUCT LINKAGE GROUP TREE                       #
#                                                                   #
#####################################################################

#####################################################################
#                                                                   #
#  FUNCTION :  construct_group_tree()                               #
#                                                                   #
#  USAGE    :  initialize group tree using BFS algorithm            #
#                                                                   #
#  ARGS     :  %Tree       ----  group tree hash                    #
#              %all_loci   ----  all loci set                       #
#              %mLOD       ----  modified LOD score hash            #
#              $start      ----  start lod                          #
#              $end        ----  end lod                            #
#              $stepSize   ----  step size of lod                   #
#                                                                   #
# RETURN    :  NONE                                                 #
#####################################################################

print $log "--- construct linkage group tree ---\n";

if (-f "$dOut/$fKey.gTree.hash" &&  $redo==0){
	
	print "$dOut/$fKey.gTree.hash already exists, retrieve group tree from file\n";

	$Tree = retrieve("$dOut/$fKey.gTree.hash");
	$succeed_flag = 1;
}else{

	construct_group_tree($Tree,\%all_loci,\%mLOD,$start,$end,$stepSize);
	store  $Tree,"$dOut/$fKey.gTree.hash";
}
print $log "\n";
out_sub_tree($Tree,"$dOut/$fKey.gTree.hash.dumper");

#####################################################################
#                                                                   #
#           SIMPLIFY GROUP TREE : DELETE GROUP FRAGMENT             #
#                                                                   #
#####################################################################

#####################################################################
#                                                                   #
#  FUNCTION :  delete_group_fragment()                              #
#                                                                   #
#  USAGE    :  delete group node whose size is less than minGROUP   #
#                                                                   #
#  ARGS     :  %Tree       ----  group tree hash                    #
#              $minGroup   ----  minimum group size                 #
#              $end        ----  end lod                            #
#                                                                   #
# RETURN    :  NONE                                                 #
#####################################################################

print $log "--- simplify group tree : delete group fragement ---\n";
if (-f "$dOut/$fKey.gTree.delFragment.hash" && $succeed_flag == 1 && $redo==0) {
	 print "$dOut/$fKey.gTree.delFragment.hash already exists, retrieve simplified group tree from file\n";

	$Tree = retrieve("$dOut/$fKey.gTree.delFragment.hash");
	$succeed_flag = 1;

}else{

	delete_group_fragment($Tree,$minGroup,$end);
	store  $Tree,"$dOut/$fKey.gTree.delFragment.hash";
}
print $log "\n";
out_sub_tree($Tree,"$dOut/$fKey.gTree.delFragment.hash.dumper");


#####################################################################
#                                                                   #
#           SIMPLIFY GROUP TREE : MERGE NODE && MERGE PATH          #
#                                                                   #
#####################################################################

#####################################################################
#                                                                   #
#  FUNCTION :  simplify_group_tree()                                #
#                                                                   #
#  USAGE    :  replace single path with root node of sub_tree       #
#                                                                   #
#  ARGS     :  %Tree       ----  group tree hash                    #
#              $end        ----  end lod                            #
#                                                                   #
# RETURN    :  NONE                                                 #
#####################################################################

print $log "--- simplify group tree : merge node && merge path ---\n";
print "--- simplify group tree : merge node && merge path ---\n";

if (-f "$dOut/$fKey.gTree.delFragment.mergeNode.hash" && $succeed_flag == 1 && $redo==0){

	print "$dOut/$fKey.gTree.delFragment.mergeNode.hash already exists, retrieve simplified group tree from file\n";

	$Tree = retrieve("$dOut/$fKey.gTree.delFragment.mergeNode.hash");
	$succeed_flag = 1;

}else{
	simplify_group_tree($Tree,$end);

	store  $Tree,"$dOut/$fKey.gTree.delFragment.mergeNode.hash";
}
print $log "\n";

out_sub_tree($Tree,"$dOut/$fKey.gTree.delFragment.mergeNode.hash.dumper");

#####################################################################
#                                                                   #
#      SELECT SUB_TREE && GET FEASIBLE SELECTION OF EACH TREE       #
#                                                                   #
#####################################################################

#####################################################################
#                                                                   #
#  FUNCTION :  select_sub_tree()                                    #
#                                                                   #
#  USAGE    :  traverse simplified tree��get a forest in which      #
#            trees has appropriate root node, calculate possible    #
#         scheme of linakge group selection for each tree in forest #
#                                                                   #
#  ARGS     :  %Tree           ----  group tree hash                #
#              $minGroup       ----  minimum group size             #
#              $maxGroup       ----  maximum group size             #
#              $end            ----  end lod                        #
#              @sub_solution   ----  all possible selection of each #
#                                    sub_tree                       #
#                                                                   #
# RETURN    :  $tree_sum       ----  number of trees in forest      #
#####################################################################

print $log "--- select sub_tree && get feasible selection of each tree ---\n";
print  "--- select sub_tree && get feasible selection of each tree ---\n";

#print Dumper \%{$Tree};

my @sub_solution = ();
my $tree_sum = select_sub_tree($Tree,$minGroup,$maxGroup,$end,\@sub_solution);

store \@sub_solution, "$dOut/$fKey.subTree.array";

#open DUMP,">$dOut/$fKey.subTree.array.dumper";
#print DUMP Dumper \%{$Tree};
#close DUMP;

if ($tree_sum > $nChro) {

	print "Trees selected exceed number of chromosomes,the reason may be: relaxed group size interval\n";
	print "Narrow the group size interval and try again.Good luck!!!\n";

	print $log "Trees selected exceed number of chromosomes,the reason may be: relaxed group size interval\n";
	print $log "Narrow the group size interval and try again.Good luck!!!\n";

	exit;
}

print $log "\n";

#####################################################################
#                                                                   #
#             PICK FEASIBLE COMBINATION OF GROUP SELECTION          #
#                                                                   #
#####################################################################

print $log "--- pick feasible combination of group selection ---\n";
print  "--- pick feasible combination of group selection ---\n";

my @feasible_combination = ();
if (!defined $sub_solution[0]||scalar @{$sub_solution[0]} == 0) {
	print $log "Trees selected result is 0\n";
	
	exit;
}
for (my $i=0;$i<@{$sub_solution[0]} ;$i++) {
	my @stack = () ;
	pick_combination (0,$i,0,\@sub_solution,\@stack,\@feasible_combination);
}

print $log "\n";

#####################################################################
#                                                                   #
#                       FIND OPTIMAL SOLUTION                       #
#                                                                   #
#####################################################################

print $log "--- find optimal solution ---\n";
print  "--- find optimal solution ---\n";

my $optimal_solution = ();
unless (@feasible_combination) {

	print "can't find solution for linkage group determination,please adjust parameters and try again.Good luck!\n";
	print $log "can't find solution for linkage group determination,please adjust parameters and try again.Good luck!\n";

	exit;
}else{

	find_optimal_solution(\@sub_solution,\@feasible_combination,\$optimal_solution);
}
print "===We are done!!===\n";##<
print $log "\n";

#####################################################################
#                                                                   #
#                             OUTPUT                                #
#                                                                   #
#####################################################################

print $log "--- output ---\n";

output_linkage_group($optimal_solution,"$dOut/$fKey");

print $log "\n";


`touch $dOut/determine_linkage_group.check`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close ($log) ;
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub out_sub_tree{
	my($tree,$file)=@_;
	my @queue=();
	push @queue,$tree->{'root'};
	open Out,">$file";
	for(my $i=0;$i<@queue;$i++) {
		next if ($queue[$i]->{'nloc'} < $minGroup);
		print Out join("\n","\t" x $queue[$i]->{'lod'}.$queue[$i]->{'lod'}. "/".$queue[$i]->{'nloc'}),"\n";
		print Out join("\n","\t" x $queue[$i]->{'lod'}.join(",",@{$queue[$i]->{'locus'}})),"\n";
		next if (!defined $queue[$i]->{'child'});
		if (scalar @{$queue[$i] ->{'child'}}!=0) {
			@queue = (@queue[0..$i],@{$queue[$i]->{'child'}},@queue[$i+1..$#queue]);
		}
	}
	close Out;
}

sub construct_group_tree {#
	my ($tree,$loci,$mLOD,$b_lod,$e_lod,$s_lod) = @_;

	my @queue = ();  ## Group node queue
	my $lod;         ## Lod level of dividing group

	## init root node 

	$tree->{'root'} = init_node([keys %{$loci}],$b_lod-$s_lod);
	push @queue,$tree->{'root'};

	## use Breadth-First-Search algorithm to construct group tree

	for(@queue){  

		$lod = $_->{'lod'} + $s_lod ;
		last if ($lod > $e_lod);

		my $child = divide_group($_,$mLOD,$lod);
		push @queue,@{$child};
		push @{$_->{'child'}},@{$child};
	}
}

sub delete_group_fragment {# delete group node whose size is smaller than the given minimum group size threshold
	my ($tree,$minG,$e_lod) = @_;

	my @queue = ();
	push @queue,$tree->{'root'};

	for(@queue){

		next if (not defined $_->{'lod'} or not exists $_->{'child'});  ## node has been removed
		last if ($_->{'lod'} > $e_lod);

		for (my $i=0;$i<@{$_->{'child'}} ;$i++) {

			if ($_->{'child'}->[$i]->{'nloc'} < $minG) {

				splice (@{$_->{'child'}},$i,1);
				$i--;
			}
		}
		push @queue,@{$_->{'child'}} if (@{$_->{'child'}});
	}

}

sub simplify_group_tree { ##  
	
	my ($tree,$e_lod) = @_;

	my @queue = ();
	push @queue,$tree->{'root'};

	for(@queue){

		last if ($_->{'lod'} > $e_lod);
		next if (not exists $_->{'child'});

		while (1) {
			
			last if (not exists $_->{'child'} or @{$_->{'child'}} > 1);

			if (exists $_->{'child'}->[0]->{'child'}) {
				$_->{'child'} = [@{$_->{'child'}->[0]->{'child'}}];
			}else{
				last;
			}
		}

		my @all_path = ();
		my @path_stack = ();

		get_all_path($_,\@path_stack,\@all_path);

		if (scalar @all_path == 1) {

			$_->{'child'} = [];
			next;
		}
		push @queue,@{$_->{'child'}};
	}
}
sub select_sub_tree {#
	my ($tree,$minG,$maxG,$e_lod,$solution) = @_;

	my $tree_num = 0;

	my @queue = ();
	push @queue,$tree->{'root'};

	for(@queue){
        next if (not exists $_->{'child'});
		last if ($_->{'lod'} > $e_lod);
		if ($_->{'nloc'} >= $minG && $_->{'nloc'} <= $maxG) { ## the current tree node satisfies the group size restriction

			$tree_num++;
			push @{$solution},get_feasible_selection($_);
			next;

		}elsif(scalar @{$_->{'child'}} != 0){ ## current node size is large than maximum group size allowed, and has child  
			push @queue,@{$_->{'child'}};
		
		}else{ 
			## current node size is large than maximum group size allowed, but has no child node
			#  this suggests that the max lod is too small to divide this group node into small children groups
			#  the program spawn a new process to perform linkage grouping in a higher lod level   
			
			#%mLOD = ();       ## hash to store pair wise data
			#$Tree = {};        ## group tree hash
			#%all_loci = ();   ## all loci set
			#_release_hash_resource(\%mLOD, \%{$tree}, \%all_loci);  

##>			print "current max lod is too small, the program restart and perform linkage grouping again\n";
##>			print $log "current max lod is too small, the program restart and perform linkage grouping again\n";

			print "current max lod is too small\n"; ##<
			print $log "current max lod is too small\n"; ##<
##>			$end += 5;  ## increase lod

##>			my $job = "perl $0 -i $fIn -k $fKey -d $dOut -n $nChro -b $start -e $end -s $stepSize -minGroup $minGroup -maxGroup $maxGroup -mode $mode -redo 1";
##>			$job .= " -g $genotype" if (defined $genotype);	##<
			
##>			print $log "$job\n";
##>			print "$job\n";
			
##>			my $log_file = _get_log_file();
			
##>			print "==WTF???=$end=\n";##<
##>			`$job >>$log_file &`;##<
##>			print "==end of WTF??=$end=\n";##<
			exit;
		}
	}
	return $tree_num;
}

## add by macx 2013-12-06
sub _get_log_file {
	return "$dOut/$fKey.linkageGroupping.log";
} 

sub _release_hash_resource{
	$_ = {} for (@_);
}

sub print_solution_matrix {#
	my ($arr) = @_;

	for (my $i=0;$i<@{$arr} ;$i++) {
		print join("\t",map {scalar @{$_}} @{$arr->[$i]}),"\n";
	}
}

sub output_linkage_group {#
	my ($final_solution,$key) = @_;		

	foreach my $nr (keys %{$final_solution}) {
			
		open (OUT,">$key.lg") or die $!;

		my $group_no = 0;

		foreach my $group (@{$final_solution->{$nr}}) {

			print OUT ">LG",++$group_no,"\tlod = ",$group->{'lod'},"\t nloc = ",$group->{'nloc'},"\n";

			print OUT join("\t",@{$group->{'locus'}}),"\n";
		}

		close (OUT) ;
	}
		
}
sub find_optimal_solution {#
	
	my ($sub_solution,$solution,$optimal) = @_;

	my %efficacy = ();
	my %lgVariance = ();

	for (my $i=0;$i<@{$solution} ;$i++) {

		my $sum = 0;
		my $sum_X2 = 0;
		my @temp = ();

		for (my $j=0;$j<@{$solution->[$i]} ;$j++) {
			next if (@{$sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]} == 0) ;

			push @temp,@{$sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]};
			for (my $k=0;$k<@{$sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]} ;$k++) {
				$sum += $sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]->[$k]->{'nloc'};
				$sum_X2 += $sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]->[$k]->{'nloc'}**2;
			}
		}

		print "sum = $sum \t plan id = $i\n";
		my $var = ($sum_X2 / @{$solution->[$i]}) - ($sum / @{$solution->[$i]})**2;
		$efficacy{$sum}{$i} = \@temp;
		$lgVariance{$var}{$i} = \@temp;
	}
##>	print "====WHY===NOT===HERE????\n";
	mkdir "$dOut/all_scheme" unless (-d "$dOut/all_scheme") ;

	foreach my $n (keys %efficacy) {
		output_linkage_group(\%{$efficacy{$n}},"$dOut/all_scheme/$fKey.$n");
	}

	if ($mode == 1) {
		my ($max) = sort {$b <=> $a} keys %efficacy;
		$$optimal = {%{$efficacy{$max}}};
	}elsif($mode == 2){
		my ($best) = sort {$a <=> $b} keys %lgVariance; 
		$$optimal = {%{$lgVariance{$best}}};
	}else{
		die "wrong mode!!!";
	}
}

sub find_optimal_solution_old {#
	my ($sub_solution,$solution,$optimal) = @_;

	my %efficacy = ();

	for (my $i=0;$i<@{$solution} ;$i++) {

		my $sum = 0;
		my @temp = ();

		for (my $j=0;$j<@{$solution->[$i]} ;$j++) {

			next if (@{$sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]} == 0) ;

			push @temp,@{$sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]};

			for (my $k=0;$k<@{$sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]} ;$k++) {
				$sum += $sub_solution->[$solution->[$i]->[$j]->[0]][$solution->[$i]->[$j]->[1]]->[$k]->{'nloc'};
			}
		}

#		print "sum = $sum \t plan id = $i\n";

		$efficacy{$sum}{$i} = \@temp;
	}

	mkdir "$dOut/all_scheme" unless (-d "$dOut/all_scheme") ;

	foreach my $n (keys %efficacy) {

		output_linkage_group(\%{$efficacy{$n}},"$dOut/all_scheme/$fKey.$n");

	}
	
	my ($max) = sort {$b <=> $a} keys %efficacy;

	$$optimal = {%{$efficacy{$max}}};

}

sub get_all_path {
	my ($node,$path,$stack) = @_;
	
	push @{$path},$node;

	if (not exists $node->{'lod'} or $node->{'lod'} == $end or  not exists $node->{'child'} or !@{$node->{'child'}}) {

#		print "leaf node\n";

		push @{$stack},[@{$path}];
		pop @{$path};
		return;
	}else{

#		print "non_leaf node\n";
		for (my $i=0;$i<@{$node->{'child'}} ;$i++) {
			get_all_path($node->{'child'}->[$i],$path,$stack);
		}
		pop @{$path};
		return;
	}

}

sub print_path {
	my ($path_arr) = @_;

	for (my $i=0;$i<@{$path_arr} ;$i++) {
	
		print ">path $i\n";
		print join("\t",map {$_->{'index'}} @{$path_arr->[$i]}),"\n";
	}

}

sub get_feasible_selection{
	my ($node) = @_;

	my @all_path = ();
	my @path_stack = ();

	get_all_path($node,\@path_stack,\@all_path);

	my @result = ();

	for (my $i=0;$i<@{$all_path[0]} ;$i++) {

		my @stack = ();
		exhaustive_search_group(0,$i,\@all_path,\@stack,\@result);

	}

	@result = @{merge_path(\@result)};  ## uniq selection

	return \@result;

}

sub exhaustive_search_group {#
	my ($i,$j,$arr,$stack,$solution) = @_;

	push @{$stack},$arr->[$i][$j];

	if ($i == @{$arr} - 1) {

		push @{$solution},[@{merge_node(\@{$stack})}];
		pop @{$stack};
		return ; 
	}else{

		for (my $k=0;$k<@{$arr->[$i+1]} ;$k++) {
			exhaustive_search_group($i+1,$k,$arr,$stack,$solution);
		}
		pop @{$stack};
		return ; 
	}
}

sub merge_node {#
	my ($arr) = @_;

	my %index = map {$_,$arr->[$_]} 0..@{$arr}-1;
	for (my $i=0;$i<@{$arr}-1 ;$i++) {
		for (my $j=$i+1;$j<@{$arr} ;$j++) {

			if (exists $index{$i} && exists $index{$j}) {
				if (isImmediateFamily([$index{$i}],$index{$j})) {
	
					if ($index{$i}->{'nloc'} > $index{$j}->{'nloc'}) {

						delete $index{$j};
					}else{

						delete $index{$i};
					}
				}
			}
		}
	}

	return [values %index];
}

sub merge_path {
	my ($arr) = @_;

	my %index = ();

	for (my $i=0;$i<@{$arr} ;$i++) {

		$index{join("|",sort { $a cmp $b } map {$_->{'index'}} @{$arr->[$i]})} = [@{$arr->[$i]}]; 
	}

	return [values %index];
}

sub isImmediateFamily {#
	my ($arr,$node) = @_;

	return 0 if (scalar @{$arr} == 0) ;

	my $isImmediateFamily = 0;

	for(@{$arr}){

		if (@{intersection([@{$_->{'locus'}}],[@{$node->{'locus'}}])}){
			 $isImmediateFamily = 1;

			 last;
		}
	}

	return $isImmediateFamily;
}

sub intersection {#
	my ($A,$B)=@_;
	my %uniqA=map {$_,1} @{$A};
	my %uniqB=map {$_,1} @{$B};
	my %merge=();
	my %overlap=();
	foreach  (keys %uniqA,keys %uniqB) {
		$merge{$_}++ && $overlap{$_}++;
	}
	my @result = keys %overlap;
	return \@result;
}

sub pick_combination {
	my ($i,$j,$sum,$arr,$stack,$answer) = @_;
	
	push @{$stack},[$i,$j];
	$sum += scalar @{$arr->[$i][$j]};

#	print "sum = ",$sum,"\n";
#	print "stack element: ",join("\t",map {scalar @{$arr->[$_->[0]][$_->[1]]}} @{$stack}),"\n";

	if ($i == @{$arr}-1) {

		if ($sum == $nChro){
#			print "solution found\n";
			push @{$answer},[@{$stack}];
			pop @{$stack};
			return ;
		}else{

			pop @{$stack};
			return ;
		}
	
	}elsif($sum > $nChro){
		pop @{$stack};
		return ;
	}else{
		for (my $k=0;$k<@{$arr->[$i+1]} ;$k++) {
			pick_combination($i+1,$k,$sum,$arr,$stack,$answer);
		}
		pop @{$stack};
		return ;
		
	}
} 

sub init_node {
	my ($arr,$lod) = @_;
	my $hash = {};
	$hash->{'lod'} = $lod;
	$hash->{'locus'} = $arr;
	$hash->{'nloc'} = scalar @{$arr};
	$hash->{'index'} = $lod."/".$hash->{'nloc'}."/".$arr->[0];
	return $hash;
}

sub divide_group {#
	my ($father,$ref_lod,$criterion) = @_;
	
	return [init_node([@{$father->{'locus'}}],$criterion)] if ($father->{'nloc'} == 1);

	my (%T,@S) = ();
	my @result =();
	

	#init
	%T = map {$_,1} @{$father->{'locus'}};
	push @S,$father->{'locus'}->[0];
	delete $T{$father->{'locus'}->[0]};

	my $last_cycle_num = -1;
	
	# interation 	
	while (1) {

		$last_cycle_num = scalar keys %T;
		my $cut_point = '';
		
		while (1) {

			foreach my $g (reverse @S) { 

				last if ($g eq $cut_point) ; # avoid repeated calculation

				foreach my $ug (keys %T) {
					if (exists $ref_lod->{$ug}{$g} and $ref_lod->{$ug}{$g} >= $criterion) {

						push @S,$ug;
						delete $T{$ug};
					
					}

					$cut_point = $g;
				}

			}

			last if (scalar keys %T == $last_cycle_num or scalar keys %T == 0) ;
			$last_cycle_num = scalar keys %T;

		}

		push @result,init_node([@S],$criterion);  ## a group node born

		last if (scalar keys %T == 0) ;

		my ($next_start) = keys %T;
		@S = ($next_start);	
		delete $T{$next_start};
	}

	return \@result;

}

sub max {
	
	my $max = shift;

	for(@_){

		$max = $_ if ($_ > $max);
	}

	return $max;
}

sub print_log_head {#
	my $fh=shift;
	my $user=`whoami`;chomp $user;
	my $program="$Bin/$Script";
	my $time=GetTime();
	my $currentPath=`pwd`;chomp $currentPath;
	my $job="perl $program"."\t".join("\t",@Original_ARGV)."\n";
	print $fh "==========================================\n";
	print $fh join("\n",(";user: $user",";workDir : $currentPath",";job: $job")),"\n";
	print $fh "==========================================\n";
}


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:

Usage:
  Options:
  -i		<file>	Input file , xxx.mLOD, forced
  -k		<str>	Key of output file,forced
  -d		<str>	Directory where output file produced,optional,default [./]

  -g		<file>	genotyping file only grouping these Markers


  -n		<int>	Species\' chromosome number,forced

  -b		<int>	Start lod,optional,default [3]
  -e		<int>	End lod,optional,default [20]
  -s		<int>	Stepsize of lod,optional,default [1]
  
  -minGroup	<int>	Minimum threshold of group size,optional,default [20]
  -maxGroup	<int>	Maximum threshold of group size,optional,default [500]
  
  -redo			Redo linkage grouping ignoring already existed file
  -dumper		Print tree hashes to dumper
  -h		Help

USAGE
	print $usage;
	exit;
}
