#!/usr/bin/env perl
# ---------------------------------------------------------
# 1. CRITICAL: SHIELD THE STDOUT CHANNEL
# ---------------------------------------------------------
# MCP uses STDOUT for JSON-RPC. We must move all Perl/Python/System 
# chatter to STDERR immediately to prevent "invalid character" errors.
my $MCP_OUT;
BEGIN {
    # Prefer vendored helper scripts and modules shipped with this repo, then
    # fall back to the user's historical global Linux_codes_SAM locations.
    require File::Basename;
    my $is_arch_dir = sub {
        my ($dir) = @_;
        return 0 unless -d $dir;
        my $name = File::Basename::basename($dir);
        return 1 if $name =~ /(?:-thread-multi|linux|gnu|darwin|MSWin32|cygwin|^x86_64|^aarch64|^arm64|^i[3-6]86)/i;
        return 1 if -d "$dir/auto";
        return 0;
    };
    my $self_dir = __FILE__;
    $self_dir =~ s{\\}{/}g;
    $self_dir =~ s{/[^/]+$}{};
    my $local_deps = ($self_dir ? "$self_dir/DiffGWASDeps" : "DiffGWASDeps");
    my $local_bin = ($self_dir ? "$self_dir/local/bin" : "local/bin");
    my $vendor_perl5 = ($self_dir ? "$self_dir/vendor/perl5" : "vendor/perl5");
    my $local_venv_bin = ($self_dir ? "$self_dir/.venv-pipeline/bin" : ".venv-pipeline/bin");
    my $local_venv_scripts = ($self_dir ? "$self_dir/.venv-pipeline/Scripts" : ".venv-pipeline/Scripts");
    my $local_python_record = ($self_dir ? "$self_dir/.venv-pipeline/.python-bin" : ".venv-pipeline/.python-bin");
    my $local_perl5 = ($self_dir ? "$self_dir/local/perl5/lib/perl5" : "local/perl5/lib/perl5");
    my @local_perl5_arch = grep { $is_arch_dir->($_) } glob("${local_perl5}/*");
    $ENV{PATH} = join(':', grep { defined && length } $local_venv_bin, $local_venv_scripts, $local_bin, $local_deps, ($ENV{PATH} // ''), '/mnt/g/NGS_lib/Linux_codes_SAM', '//rs1.stjude.org/clusterhome/zcheng/NGS_lib/Linux_codes_SAM');
    $ENV{PERL5LIB} = join(':', grep { defined && length } $vendor_perl5, $local_perl5, @local_perl5_arch, $local_deps, ($ENV{PERL5LIB} // ''), '/mnt/g/NGS_lib/Linux_codes_SAM', '//rs1.stjude.org/clusterhome/zcheng/NGS_lib/Linux_codes_SAM');
    if (!$ENV{PIPELINE_PYTHON_BIN} && -f $local_python_record) {
        if (open(my $pfh, '<', $local_python_record)) {
            my $line = <$pfh>;
            close $pfh;
            if (defined $line) {
                chomp $line;
                $ENV{PIPELINE_PYTHON_BIN} = $line if length $line;
            }
        }
    }
    open($MCP_OUT, ">&STDOUT") or die "Could not duplicate STDOUT: $!";
    open(STDOUT, ">&STDERR") or die "Could not redirect STDOUT: $!";
    select(STDERR); $| = 1; # Unbuffer STDERR for real-time logging
    $ENV{MOJO_LOG_LEVEL} //= 'error'; # Mute Mojolicious info logs
    select(STDERR); $| = 1;
}
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/DiffGWASDeps";
use SAS_ODA_Runner; # Ensure this module is in your Perl library path
use Mojolicious::Lite -signatures;
use MCP::Server;
use POSIX qw(strftime);
use Cwd qw(getcwd);
use HTTP::Tiny; # Core module - always available
use File::Temp qw(tempfile);
use File::Basename;
use File::Path qw(make_path);
use File::Which qw(which); # Ensure this is installed, or use the backtick version below
use Mojo::Util   qw(md5_sum);
use JSON::MaybeXS qw(encode_json);

sub shell_quote_single {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/'/'"'"'/g;
    return "'$text'";
}

sub powershell_quote_single {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/'/''/g;
    return "'$text'";
}

sub cygpath_to_windows {
    my ($path) = @_;
    return $path unless defined $path && length $path;
    return $path if $path =~ m{^[A-Za-z]:[\\/]} || $path =~ m{^\\\\};
    my $quoted = shell_quote_single($path);
    my $win = `cygpath -w $quoted`;
    chomp $win;
    return length($win) ? $win : $path;
}

sub cleanup_generated_tmpdir_if_empty {
    my ($path) = @_;
    return unless defined $path && length $path;

    my $dir = -d $path ? $path : dirname($path);
    return unless defined $dir && length $dir && -d $dir;

    my $base = basename($dir);
    return unless defined $base && $base =~ /^tmp\d+$/;

    opendir(my $dh, $dir) or return;
    my @entries = grep { $_ ne '.' && $_ ne '..' } readdir($dh);
    closedir $dh;
    return if @entries;

    rmdir $dir;
}

sub cleanup_all_generated_empty_tmpdirs {
    for my $dir (glob('tmp*')) {
        next unless defined $dir && -d $dir;
        cleanup_generated_tmpdir_if_empty($dir);
    }
}

sub pid_is_running {
    my ($pid) = @_;
    return 0 unless defined $pid && $pid =~ /^\d+$/;

    my @pids = split(" ", `perl -S GetPIDs.pl`);
    chomp @pids;
    return 1 if grep { $_ == $pid } @pids;

    if ($^O =~ /^(?:cygwin|MSWin32)$/i) {
        my $task = `tasklist /FI "PID eq $pid" /FO CSV /NH 2>NUL`;
        return 0 unless defined $task && length $task;
        return 0 if $task =~ /No tasks are running/i;
        return 1 if $task =~ /^"[^"]+"/;
    }

    return 0;
}

if (@ARGV==0) {
    print <<USAGE;
Usage: server.pl daemon -m production -l http://127.0.0.1:8080 
This script will run a Perl MCP server accessible by AI agent locally;
It is intended for local editor-driven workflows such as VS Code plus Codex,
including Windows/Cygwin, macOS, and Linux shells. When Ubuntu validation
cannot be completed through Vagrant, the same workflow has also been tested
through an isolated Docker Ubuntu runtime, and the repository now also ships a
saved Dockerfile plus Singularity/Apptainer definition for containerized use.
USAGE
    exit 0;
}


my %sessions;

my $server = MCP::Server->new;

$server->tool(
    name        => 'run_perl_or_bash_cmd',
    description => 'run perl_or_bash codes or script by submitting the codes/scripts as a background job, and return the results in a text file '.
    'containing enough information of perl_or_bash log and others for debugging later. This tool is designed to avoid timeout issues when running ' .
    'long queries on platforms like Gemini. The tool checks the status of the background job and returns partial output if still running,'.
    'or the final results when complete.',
    input_schema => {
        type => 'object',
        properties => {
            perl_or_bash_codes_or_file  => { type => 'string', 
                                   description => 'The RAW, UNMODIFIED perl_or_bash code or file path. ' .
                           'IMPORTANT: You must preserve all spaces, newlines, semicolons, and special characters exactly. ' .
                           'Do not perform any text normalization, "cleaning," or whitespace removal.'},
            output_file   => { type => 'string',
                               description => 'Optional output file prefix for output files containing perl_or_bash log and other useful '.
                               'information for debugging (optional); default output_file name is "output", '. 
                               'which will be added with the appendix ".html.info.txt", all the final default output file, "./tmp*/output.html.info.txt", '. 
                               'will be put into a temporary directory named as "./tmp*", with * respresenting 6 random numbers.' },
            is_bash       => { type=>'string',
                               description => 'Optional for running the input command using bash -c, otherwise, using perl -S'},
            pid           => { type => 'integer',
                               description => 'Optional PID to check status of a previous query (when run the tool again with the same perl_or_bash_codes_or_file)' },
            tmp_perl_or_bash_file     => { type => 'string',
                               description => 'Optional internal use only: temporary perl_or_bash file path to store the perl_or_bash codes when the input perl_or_bash codes contain' . 
                               'special characters that may cause issues with passing the perl_or_bash codes as a command line argument to the python function. '. 
                               'The tool will save the perl_or_bash codes into this temporary perl_or_bash file and run the query with the temporary perl_or_bash file, and then ' . 
                               'delete the temporary perl_or_bash file after running the query.' }                   
        },
        required => ['perl_or_bash_codes_or_file']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        #1. Debug: Check if the spaces exist the moment Perl gets them
        #(Using a log file because STDOUT is for the MCP protocol)
        # open my $log, '>>', 'output.html.info_debug.log';
        # my $raw_input = $args->{perl_or_bash_codes_or_file};
        # print $log "RECEIVED: [$raw_input]\n";
        # close $log;   

        my $perl_or_bash_codes_or_file = $args->{perl_or_bash_codes_or_file};
        my $tmp_perl_or_bash_file = $args->{tmp_perl_or_bash_file}; # This is for storing the temporary perl_or_bash file path when the perl_or_bash codes contain special characters;
        my $run_it_with_bash=$args->{is_bash} // 1;
        my $script_appendix='.sh';
           $script_appendix='.pl' unless $run_it_with_bash;
        #need to double quote the perl_or_bash codes when passing to the python function, 
        #otherwise the special characters and spaces may cause issues with the parsing of the input in the python function;
        my $perl_or_bash_codes;
        if (length($perl_or_bash_codes_or_file)<50 && -f "$perl_or_bash_codes_or_file"){
            $perl_or_bash_codes=$perl_or_bash_codes_or_file; 
            }else{
            #it would be prone to error if the perl_or_bash codes contains dobule quotes;
            #To resolve the issue, we can escape the double quotes in the perl_or_bash codes before passing it to the python function;
            $perl_or_bash_codes = $perl_or_bash_codes_or_file;
            if ($perl_or_bash_codes =~ /("|')/) {#write these codes into a perl_or_bash file to run it, to avoid the issue caused by escaping the quotes in the perl_or_bash codes;
                $tmp_perl_or_bash_file = "./tmp_perl_or_bash_codes_" . time() . "$script_appendix";#Need to add ./;
                open(my $fh, '>', $tmp_perl_or_bash_file) or die "Could not create temporary bash or perl file: $!";
                print $fh $perl_or_bash_codes;
                close $fh;
                `chmod a+x $tmp_perl_or_bash_file`;
                $perl_or_bash_codes=$tmp_perl_or_bash_file;
            } else {
                $perl_or_bash_codes="$perl_or_bash_codes";
            }
        }
        
        if ($run_it_with_bash){
               if ($perl_or_bash_codes=~/^\.\/\S+\.sh/){
                $perl_or_bash_codes="bash $perl_or_bash_codes";
               }else{
                $perl_or_bash_codes="bash -c \"$perl_or_bash_codes\" ";
               }
              
            }else{
                $perl_or_bash_codes="perl -S ".'\"'.$perl_or_bash_codes.'\"';
         }
        #This output filename if provided with the previous pid, the tool will check the status of the previous query
        #and return the results in this file when the job is finished. If not provided, the tool will save the results 
        #in a default file named "tmp{timestamp}/haploreg_{query_snp}.tsv" as demonstrated in the later part of the code;
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        #Note: it is necessary to define the $out_file variable here to avoid 
        #the "Use of uninitialized value $out_file in concatenation (.) or 
        #string at line codes" error when the tool is called for the first time without the pid argument, 
        #as the $out_file variable is used in the return statement for the initial query status. 
        #Defining it here with an undefined value allows us to use it in the return statement 
        #without causing an error, and it will be properly assigned later when checking the PID or starting a new query.
        my $out_file;
        my $pid_file;
        # If checking a previous PID
        if (defined $pid_arg) {
            $pid_file = "tmp*/output.info.pid";
            my @pid_files = glob($pid_file);
            
            foreach my $pf (@pid_files) {
                if (-f $pf) {
                    open(my $pfh, '<', $pf);
                    my $stored_pid = <$pfh>;
                    close $pfh;
                    chomp $stored_pid;
                    
                    if ($stored_pid == $pid_arg) {
                        $out_file = $pf;
                        #Note: the output file is named "tmp{timestamp}/perl_or_bash_ODA_html_output_name.tsv" by default when the query is first 
                        #run without the output_file argument, and the output file path is saved in the same directory as the pid file. 
                        #So when checking the status with the pid, we can get the output file path by replacing the ".pid" extension 
                        #in the pid file name with ".txt". This way, we can ensure that we are checking the correct output file 
                        #associated with the specific PID.
                        $out_file =~ s/\.pid$/.txt/;
                        
                        my @pids = split(" ", `perl -S GetPIDs.pl`);
                        chomp @pids;
                        
                        if (grep { $_ == $pid_arg } @pids) {
                            return {
                                content => [{
                                    type => "text",
                                    text => "STATUS: RUNNING (PID $pid_arg)\n" .
                                           "Output file: $out_file\nAsk the AI agent to check status again in a moment."
                                }]
                            };
                        }
                        
                        # Job finished
                        my $content = "";
                        if (-f $out_file && open(my $fh, '<', $out_file)) {
                            local $/;
                            $content = <$fh>;
                            close $fh;
                        }
                        
                        unlink $pf;
                        cleanup_generated_tmpdir_if_empty($pf);
                        
                        return {
                            content => [{
                                type => "text",
                                text => "STATUS: COMPLETE (PID $pid_arg)\n\nperl_or_bash log for debugging saved to: $out_file\n\n" . $content
                            }]
                        };
                    }
                }
            }
            
            return {
                content => [{
                    type => "text",
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        # Start new query by creating a tmp directory and running the query in the background
        my $tmpdir = "tmp" . time();
        unless (-d $tmpdir) {
            mkdir $tmpdir or die "Failed to create tmp directory: $!";      
        }
        # Create a file named "tmp{timestamp}/haploreg_{query_snp}.tsv" for first time query, 
        #and save the results in this file when the job is finished. If the tool is called again 
        #with the same query_snp and the pid of the previous job, it will check the status of the 
        #previous job and return the results in this file when the job is finished;
        $pid_file = "$tmpdir/output.info.pid";

        #Note: the output file is named "tmp{timestamp}/output.html.info.txt" by default 
        #when the query is first run without the output_file argument, and the output file 
        #path is saved in the same directory as the pid file.
        $out_file = $output_file // "$tmpdir/output.info.txt";
        
        #Note: the temporary perl_or_bash file is named "tmp_perl_or_bash_codes_{timestamp}.perl_or_bash" when the perl_or_bash codes
        # contain special characters and need to be saved into a temporary perl_or_bash file to run the query, 
        #and the temporary perl_or_bash file is saved in the same directory as the pid file.
        $tmp_perl_or_bash_file = $args->{tmp_perl_or_bash_file} // "$tmpdir/tmp_perl_or_bash_file$script_appendix"; # This is for storing the temporary perl_or_bash file path when the perl_or_bash codes contain special characters;
        
        unless (-f $pid_file) {
            my $pid = fork();
            return { content => [{ type => "text", text => "ERROR: Could not fork." }] } unless defined $pid;

            if ($pid == 0) {
                my $cmd = "$perl_or_bash_codes >$out_file";
                print STDERR "Executing command in background: $cmd\n";
                exec "$cmd 2>&1";
                exit(1);
            }

            open(my $pfh, '>', $pid_file);
            print $pfh $pid;
            close $pfh;
            # Return initial status with PID and instructions for checking later, avoiding the endless loop of checking by the AI agent, 
            # which may cause excessive checking and potential issues with the server;
            # If the perl_or_bash_codes_or_file contains single/double quotes, the tool will save the perl_or_bash codes into a temporary perl_or_bash file;
            if (-f $tmp_perl_or_bash_file){
              return {
                content => [{
                    type => "text",
                    text => "QUERYING: perl_or_bash ODA for $perl_or_bash_codes\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"perl_or_bash_codes\": \"$tmp_perl_or_bash_file\", \"pid\": $pid}\n" . 
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"
                }]
            };
            }else{
              return {
                content => [{
                    type => "text",
                    text => "QUERYING: perl_or_bash ODA for $perl_or_bash_codes_or_file\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"perl_or_bash_codes_or_file\": \"$perl_or_bash_codes_or_file\", \"pid\": $pid}\n" .
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"}]
        };
       }
    }
  }
);

$server->tool(
    name        => 'run_sas_codes_or_script_on_local_Windows',
    description => 'run sas codes or script by submitting the codes/scripts as a background job to local Windows SAS, and return SAS log in a text file, called output.html.info.log. ' .
    'When input sas codes, the tool will automatically save these codes into a temperary sas script and run it using local SAS.' .
    'It will return the results in a text file, containing enough information of SAS log and others for debugging later. This tool is designed to avoid timeout issues when running ' .
    'long queries on platforms like Gemini. The tool checks the status of the background job and returns partial output if still running,' .
    'or the final results when complete.',
    input_schema => {
        type => 'object',
        properties => {
            sas_codes_or_file  => { type => 'string', 
                                   description => 'The RAW, UNMODIFIED SAS code or file path. ' .
                           'IMPORTANT: You must preserve all spaces, newlines, semicolons, and special characters exactly. ' .
                           'Do not perform any text normalization, "cleaning," or whitespace removal.'},
            output_file   => { type => 'string',
                               description => 'Optional output file prefix for output files containing SAS log and other useful '.
                               'information for debugging (optional); default output_file name is "output", '. 
                               'which will be added with the appendix ".html.info.txt", all the final default output file, "./tmp*/output.html.info.txt", '. 
                               'will be put into a temporary directory named as "./tmp*", with * respresenting 6 random numbers.' },                           
            pid           => { type => 'integer',
                               description => 'Optional PID to check status of a previous query (when run the tool again with the same sas_codes_or_file)' },
            tmp_sas_file     => { type => 'string',
                               description => 'Optional internal use only: temporary SAS file path to store the SAS codes when the input are sas codes. ' .
                               'The tool will save the sas codes into this temporary sas file and run the query with the temporary sas file, and then ' . 
                               'delete the temporary sas file after running the query.' }                   
        }
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        #1. Debug: Check if the spaces exist the moment Perl gets them
        #(Using a log file because STDOUT is for the MCP protocol)
        # open my $log, '>>', 'output.html.info_debug.log';
        # my $raw_input = $args->{sas_codes_or_file};
        # print $log "RECEIVED: [$raw_input]\n";
        # close $log;   

        my $sas_codes_or_file = $args->{sas_codes_or_file} // ' '; # Default to a single space if not provided, to avoid issues with empty input
        my $tmp_sas_file = $args->{tmp_sas_file}; # This is for storing the temporary sas file path when the sas codes contain special characters;
        #need to double quote the sas codes when passing to the python function, 
        #otherwise the special characters and spaces may cause issues with the parsing of the input in the python function;
        my $sas_codes="";
        if (-f "$sas_codes_or_file"){
            $sas_codes="$sas_codes_or_file"; 
            }else{
            #it would be prone to error if the sas codes contains dobule quotes;
            #To resolve the issue, we can escape the double quotes in the sas codes before passing it to the python function;
            $sas_codes = $sas_codes_or_file;
            $tmp_sas_file = "tmp_sas_codes_" . time() . ".sas";
            open(my $fh, '>', $tmp_sas_file) or die "Could not create temporary SAS file: $!";
            print $fh $sas_codes;
            close $fh;
            $sas_codes="$tmp_sas_file";
        }

        #Need to add ods html at the first line to enable the final html will be printed;
        open F,"$sas_codes" or die "failed to open the file $sas_codes: $!";
        open O,">$sas_codes.new" or die "failed to update the file $sas_codes.new: $!";
        while (my $l=<F>){
        print O "ods html path='.';\n" if $.==1 and $l !~/^ods html path\=/;
        print O $l;
        }
        close F;
        close O;
        unlink $sas_codes;
        `mv $sas_codes.new $sas_codes`;

        
        #This output filename if provided with the previous pid, the tool will check the status of the previous query
        #and return the results in this file when the job is finished. If not provided, the tool will save the results 
        #in a default file named "tmp{timestamp}/haploreg_{query_snp}.tsv" as demonstrated in the later part of the code;
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        #Note: it is necessary to define the $out_file variable here to avoid 
        #the "Use of uninitialized value $out_file in concatenation (.) or 
        #string at line codes" error when the tool is called for the first time without the pid argument, 
        #as the $out_file variable is used in the return statement for the initial query status. 
        #Defining it here with an undefined value allows us to use it in the return statement 
        #without causing an error, and it will be properly assigned later when checking the PID or starting a new query.
        my $out_file;
        my $pid_file;
        # If checking a previous PID
        if (defined $pid_arg) {
            $pid_file = "tmp*/output.html.info.pid";
            my @pid_files = glob($pid_file);
            
            foreach my $pf (@pid_files) {
                if (-f $pf) {
                    open(my $pfh, '<', $pf);
                    my $stored_pid = <$pfh>;
                    close $pfh;
                    chomp $stored_pid;
                    
                    if ($stored_pid == $pid_arg) {
                        $out_file = $pf;
                        #Note: the output file is named "tmp{timestamp}/SAS_ODA_html_output_name.tsv" by default when the query is first 
                        #run without the output_file argument, and the output file path is saved in the same directory as the pid file. 
                        #So when checking the status with the pid, we can get the output file path by replacing the ".pid" extension 
                        #in the pid file name with ".txt". This way, we can ensure that we are checking the correct output file 
                        #associated with the specific PID.
                        $out_file =~ s/\.pid$/.txt/;
                        
                        my @pids = split(" ", `perl -S GetPIDs.pl`);
                        chomp @pids;
                        
                        if (grep { $_ == $pid_arg } @pids) {
                            return {
                                content => [{
                                    type => "text",
                                    text => "STATUS: RUNNING (PID $pid_arg)\n" .
                                           "Output file: $out_file\nAsk the AI agent to check status again in a moment."
                                }]
                            };
                        }
                        
                        # Job finished
                        my $content = "";
                        if (-f $out_file && open(my $fh, '<', $out_file)) {
                            local $/;
                            $content = <$fh>;
                            close $fh;
                        }
                        
                        unlink $pf;
                        cleanup_generated_tmpdir_if_empty($pf);
                        
                        return {
                            content => [{
                                type => "text",
                                text => "STATUS: COMPLETE (PID $pid_arg)\n\nSAS log for debugging saved to: $out_file\n\n" . $content
                            }]
                        };
                    }
                }
            }
            
            return {
                content => [{
                    type => "text",
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        # Start new query by creating a tmp directory and running the query in the background
        my $tmpdir = "tmp" . time();
        unless (-d $tmpdir) {
            mkdir $tmpdir or die "Failed to create tmp directory: $!";      
        }
        # Create a file named "tmp{timestamp}/haploreg_{query_snp}.tsv" for first time query, 
        #and save the results in this file when the job is finished. If the tool is called again 
        #with the same query_snp and the pid of the previous job, it will check the status of the 
        #previous job and return the results in this file when the job is finished;
        $pid_file = "$tmpdir/output.html.info.pid";

        #Note: the output file is named "tmp{timestamp}/output.html.info.txt" by default 
        #when the query is first run without the output_file argument, and the output file 
        #path is saved in the same directory as the pid file.
        $out_file = $output_file // "$tmpdir/output.html.info.log";
        
        #Note: the temporary sas file is named "tmp_sas_codes_{timestamp}.sas" when the sas codes
        # contain special characters and need to be saved into a temporary sas file to run the query, 
        #and the temporary sas file is saved in the same directory as the pid file.
        $tmp_sas_file = $args->{tmp_sas_file} // "$tmpdir/tmp_sas_file.sas"; # This is for storing the temporary sas file path when the sas codes contain special characters;
        
        unless (-f $pid_file) {
            my $pid = fork();
            return { content => [{ type => "text", text => "ERROR: Could not fork." }] } unless defined $pid;

            if ($pid == 0) {
                my $cmd = "RunWindowsSAS.sh $sas_codes output.html.info $tmpdir";
                print STDERR "Executing command in background: $cmd\n";
                exec "$cmd 2>&1";
                exit(1);
            }

            open(my $pfh, '>', $pid_file);
            print $pfh $pid;
            close $pfh;
            # Return initial status with PID and instructions for checking later, avoiding the endless loop of checking by the AI agent, 
            # which may cause excessive checking and potential issues with the server;
            # If the sas_codes_or_file contains single/double quotes, the tool will save the sas codes into a temporary sas file;
            if (-f $tmp_sas_file){
              return {
                content => [{
                    type => "text",
                    text => "QUERYING: SAS ODA for $sas_codes\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"sas_codes\": \"$tmp_sas_file\", \"pid\": $pid}\n" . 
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"
                }]
            };
            }else{
              return {
                content => [{
                    type => "text",
                    text => "QUERYING: SAS ODA for $sas_codes_or_file\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"sas_codes_or_file\": \"$sas_codes_or_file\", \"pid\": $pid}\n" .
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"}]
        };
       }
    }
  }
);

$server->tool(
    name        => 'run_sas_codes_or_script_in_ODA',
    description => 'submit sas codes/scripts or upload file, delete file, list files into SAS ODA as a background job, and return the results in a text file. ' .
    'It will automatically upload a required file to remote SAS ODA with or without running codes or scriptif specified as input.' .
    'This means the tool can run sas codes or script with or without uploading a file or just upload a file, delete a file, or list files in a directory in ODA each time without running sas codes or script. ' .
    'Additionally, several extra arguments, including delete_file, delete_file_rgx, download_file, file_info, and dir4listing can work similarly like upload_file, meaning to work individually or combinationally with other arguments. ' .
    'Bulk file operations are now supported: upload_file, download_file, download_local_path, delete_file, delete_file_rgx, and file_info may each be passed either as one string or as an array of strings. ' .
    'Regex-based deletes now match both the basename and the resolved remote path, so patterns such as ~\\/.*\\.png or .*\\.png$ both work.' .
    'Remote home-directory paths should be passed literally, for example quoted as ~/plot.png from a shell layer, and the helper now normalizes both ~/... and absolute SAS home paths such as /home/... consistently across download, delete, and file-info operations. ' .
    'Delete requests are also verified after the helper reports success so stale path-mismatch cases are easier to catch. ' .
    'When persistent session reuse is enabled, uploads, downloads, deletes, directory listing, and SAS code submission now all share the same SAS ODA session instead of silently creating fresh saspy sessions for file operations. ' .
    'On the first SAS ODA run, the vendored helper now validates the saved or newly supplied SAS ODA account/password with proc setinit;run; and stores successful credentials in the SASPy authinfo file so later runs do not need to prompt again. ' .
    'The dir4listing option can be run independently when it is necessary just to check whether specific macros or files are available in a remote SAS ODA directory (default to be the HOME directory represented by ~) before running any SAS codes or file!' .
    'It will return the results in a text file, containing enough information of SAS log and others for debugging later. This tool is designed to avoid timeout issues when running ' .
    'long queries on platforms like Gemini. The tool checks the status of the background job and returns partial output if still running,' .
    'or the final results when complete.',
    input_schema => {
        type => 'object',
        properties => {
            sas_codes_or_file  => { type => 'string', 
                                   description => 'The RAW, UNMODIFIED SAS code or file path. ' .
                           'IMPORTANT: You must preserve all spaces, newlines, semicolons, and special characters exactly. ' .
                           'Do not perform any text normalization, "cleaning," or whitespace removal.'},
            output_file   => { type => 'string',
                               description => 'Optional output file prefix for output files containing SAS log and other useful '.
                               'information for debugging (optional); default output_file name is "output", '. 
                               'which will be added with the appendix ".html.info.txt", all the final default output file, "./tmp*/output.html.info.txt", '. 
                               'will be put into a temporary directory named as "./tmp*", with * respresenting 6 random numbers.' },
            upload_file   => {
                               oneOf => [
                                 { type => 'string' },
                                 { type => 'array', items => { type => 'string' } }
                               ],
                               description => 'Optional local file path or list of file paths to upload to remote SAS ODA. Each uploaded file keeps its basename in SAS ODA HOME unless the downstream wrapper renames it first.'
                             },
            delete_file   => {
                               oneOf => [
                                 { type => 'string' },
                                 { type => 'array', items => { type => 'string' } }
                               ],
                               description => 'Optional remote file path or list of remote file paths to delete from SAS ODA. Quote ~/... paths when calling through shell layers; the helper now normalizes both ~/... and absolute /home/... SAS home paths before delete and verifies that the target no longer resolves afterward.'
                             },
            delete_file_rgx => {
                               oneOf => [
                                 { type => 'string' },
                                 { type => 'array', items => { type => 'string' } }
                               ],
                               description => 'Optional regex pattern or list of regex patterns used to delete multiple remote files after listing delete_dir. Regexes match both bare filenames and resolved remote paths such as ~/plot.png.'
                             },
            delete_dir => { type => 'string',
                               description => 'Optional remote directory to scan before applying delete_file_rgx. Default: ~' },
            download_file => {
                               oneOf => [
                                 { type => 'string' },
                                 { type => 'array', items => { type => 'string' } }
                               ],
                               description => 'Optional remote file path or list of remote file paths to download from SAS ODA to the local computer. Quote ~/... paths when calling through shell layers; the helper now normalizes both ~/... and absolute /home/... SAS home paths before download.'
                             },
            download_local_path => {
                               oneOf => [
                                 { type => 'string' },
                                 { type => 'array', items => { type => 'string' } }
                               ],
                               description => 'Optional explicit local destination path or list of destination paths paired positionally with download_file.'
                             },
            file_info => {
                               oneOf => [
                                 { type => 'string' },
                                 { type => 'array', items => { type => 'string' } }
                               ],
                               description => 'Optional remote file path or list of paths for existence and size lookup. Quote ~/... paths when calling through shell layers; the helper now normalizes both ~/... and absolute /home/... SAS home paths before lookup.'
                             },
            dir4listing   => {type => 'string',
                              description => 'Optional for listing files in a remote SAS ODA directory. For the HOME directory represented by ~, pass the literal path quoted through the shell layer, for example ~/ enclosed in single quotes.'},                             
            persistent    => { type => 'string',
                               description => 'Optional truthy flag to reuse a persistent SAS ODA session across multiple tool calls. Accepts values like 1, true, yes. Default: true.' },
            session_id    => { type => 'string',
                               description => 'Optional persistent SAS ODA session id. Default: mysession.' },
            sas_oda_account => { type => 'string',
                               description => 'Optional SAS ODA account/email for noninteractive first-run login bootstrap.' },
            sas_oda_password => { type => 'string',
                               description => 'Optional SAS ODA password for noninteractive first-run login bootstrap.' },
            prompt_sas_oda_auth => { type => 'string',
                               description => 'Optional truthy flag to force a SAS ODA credential refresh even when a saved authinfo entry already exists.' },
            pid           => { type => 'integer',
                               description => 'Optional PID to check status of a previous query (when run the tool again with the same sas_codes_or_file)' },
            tmp_sas_file     => { type => 'string',
                               description => 'Optional internal use only: temporary SAS file path to store the SAS codes when the input sas codes contain' . 
                               'special characters that may cause issues with passing the sas codes as a command line argument to the python function. '. 
                               'The tool will save the sas codes into this temporary sas file and run the query with the temporary sas file, and then ' . 
                               'delete the temporary sas file after running the query.' }                   
        }
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        #1. Debug: Check if the spaces exist the moment Perl gets them
        #(Using a log file because STDOUT is for the MCP protocol)
        # open my $log, '>>', 'output.html.info_debug.log';
        # my $raw_input = $args->{sas_codes_or_file};
        # print $log "RECEIVED: [$raw_input]\n";
        # close $log;   

        my $sas_codes_or_file = $args->{sas_codes_or_file} // ' '; # Default to a single space if not provided, to avoid issues with empty input
        my $tmp_sas_file = $args->{tmp_sas_file}; # This is for storing the temporary sas file path when the sas codes contain special characters;
        my @uploaders = ref($args->{upload_file}) eq 'ARRAY'
          ? @{$args->{upload_file}}
          : (defined $args->{upload_file} ? ($args->{upload_file}) : ());
        my @downloaders = ref($args->{download_file}) eq 'ARRAY'
          ? @{$args->{download_file}}
          : (defined $args->{download_file} ? ($args->{download_file}) : ());
        my @download_local_paths = ref($args->{download_local_path}) eq 'ARRAY'
          ? @{$args->{download_local_path}}
          : (defined $args->{download_local_path} ? ($args->{download_local_path}) : ());
        my @deleters = ref($args->{delete_file}) eq 'ARRAY'
          ? @{$args->{delete_file}}
          : (defined $args->{delete_file} ? ($args->{delete_file}) : ());
        my @delete_rgxs = ref($args->{delete_file_rgx}) eq 'ARRAY'
          ? @{$args->{delete_file_rgx}}
          : (defined $args->{delete_file_rgx} ? ($args->{delete_file_rgx}) : ());
        my @file_infos = ref($args->{file_info}) eq 'ARRAY'
          ? @{$args->{file_info}}
          : (defined $args->{file_info} ? ($args->{file_info}) : ());
        my $ODADir=$args->{dir4listing};
        my $delete_dir = $args->{delete_dir};
        my $persistent_arg = $args->{persistent};
        my $session_id_arg = $args->{session_id} // 'mysession';
        my $sas_oda_account = $args->{sas_oda_account} // '';
        my $sas_oda_password = $args->{sas_oda_password} // '';
        my $prompt_sas_oda_auth = $args->{prompt_sas_oda_auth} // '';
        my $use_persistent = 1;
        if (defined $persistent_arg && $persistent_arg !~ /^(?:1|true|yes|y)$/i) {
            $use_persistent = 0;
        }
        #need to double quote the sas codes when passing to the python function, 
        #otherwise the special characters and spaces may cause issues with the parsing of the input in the python function;
        my $sas_codes="";
        if (-f "$sas_codes_or_file"){

            $sas_codes="--file \'$sas_codes_or_file\'"; 
            }else{
            #it would be prone to error if the sas codes contains dobule quotes;
            #To resolve the issue, we can escape the double quotes in the sas codes before passing it to the python function;
            $sas_codes = $sas_codes_or_file;
            if ($sas_codes =~ /("|')/) {#write these codes into a sas file to run it, to avoid the issue caused by escaping the quotes in the sas codes;
                $tmp_sas_file = "tmp_sas_codes_" . time() . ".sas";
                open(my $fh, '>', $tmp_sas_file) or die "Could not create temporary SAS file: $!";
                print $fh $sas_codes;
                close $fh;
                $sas_codes="--file \'$tmp_sas_file\'";
            } else {
                $sas_codes="--code \'$sas_codes\'";
            }
        }

        $sas_codes .= join('', map { " --upload-file " . shell_quote_single($_) } @uploaders) if @uploaders;
        $sas_codes .= join('', map { " --download-file " . shell_quote_single($_) } @downloaders) if @downloaders;
        $sas_codes .= join('', map { " --download-local-path " . shell_quote_single($_) } @download_local_paths) if @download_local_paths;
        $sas_codes .= join('', map { " --delete-file " . shell_quote_single($_) } @deleters) if @deleters;
        $sas_codes .= join('', map { " --delete-file-rgx " . shell_quote_single($_) } @delete_rgxs) if @delete_rgxs;
        $sas_codes .= " --delete-dir " . shell_quote_single($delete_dir) if defined $delete_dir && length $delete_dir;
        $sas_codes .= join('', map { " --file-info " . shell_quote_single($_) } @file_infos) if @file_infos;
        $sas_codes=$sas_codes. " --dir4listing \'$ODADir\' " if $ODADir; #If there is a ODADir as input, add this argument to listing files in remote ODA directory;
        
        #This output filename if provided with the previous pid, the tool will check the status of the previous query
        #and return the results in this file when the job is finished. If not provided, the tool will save the results 
        #in a default file named "tmp{timestamp}/haploreg_{query_snp}.tsv" as demonstrated in the later part of the code;
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        #Note: it is necessary to define the $out_file variable here to avoid 
        #the "Use of uninitialized value $out_file in concatenation (.) or 
        #string at line codes" error when the tool is called for the first time without the pid argument, 
        #as the $out_file variable is used in the return statement for the initial query status. 
        #Defining it here with an undefined value allows us to use it in the return statement 
        #without causing an error, and it will be properly assigned later when checking the PID or starting a new query.
        my $out_file;
        my $pid_file;
        # If checking a previous PID
        if (defined $pid_arg) {
            $pid_file = "tmp*/output.html.info.pid";
            my @pid_files = glob($pid_file);
            
            foreach my $pf (@pid_files) {
                if (-f $pf) {
                    open(my $pfh, '<', $pf);
                    my $stored_pid = <$pfh>;
                    close $pfh;
                    chomp $stored_pid;
                    
                    if ($stored_pid == $pid_arg) {
                        $out_file = $pf;
                        #Note: the output file is named "tmp{timestamp}/SAS_ODA_html_output_name.tsv" by default when the query is first 
                        #run without the output_file argument, and the output file path is saved in the same directory as the pid file. 
                        #So when checking the status with the pid, we can get the output file path by replacing the ".pid" extension 
                        #in the pid file name with ".txt". This way, we can ensure that we are checking the correct output file 
                        #associated with the specific PID.
                        $out_file =~ s/\.pid$/.txt/;
                        
                        my @pids = split(" ", `perl -S GetPIDs.pl`);
                        chomp @pids;
                        
                        if (grep { $_ == $pid_arg } @pids) {
                            return {
                                content => [{
                                    type => "text",
                                    text => "STATUS: RUNNING (PID $pid_arg)\n" .
                                           "Output file: $out_file\nAsk the AI agent to check status again in a moment."
                                }]
                            };
                        }
                        
                        # Job finished
                        my $content = "";
                        if (-f $out_file && open(my $fh, '<', $out_file)) {
                            local $/;
                            $content = <$fh>;
                            close $fh;
                        }
                        
                        unlink $pf;
                        cleanup_generated_tmpdir_if_empty($pf);
                        
                        return {
                            content => [{
                                type => "text",
                                text => "STATUS: COMPLETE (PID $pid_arg)\n\nSAS log for debugging saved to: $out_file\n\n" . $content
                            }]
                        };
                    }
                }
            }
            
            return {
                content => [{
                    type => "text",
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        # Start new query by creating a tmp directory and running the query in the background
        my $tmpdir = "tmp" . time();
        unless (-d $tmpdir) {
            mkdir $tmpdir or die "Failed to create tmp directory: $!";      
        }
        # Create a file named "tmp{timestamp}/haploreg_{query_snp}.tsv" for first time query, 
        #and save the results in this file when the job is finished. If the tool is called again 
        #with the same query_snp and the pid of the previous job, it will check the status of the 
        #previous job and return the results in this file when the job is finished;
        $pid_file = "$tmpdir/output.html.info.pid";

        #Note: the output file is named "tmp{timestamp}/output.html.info.txt" by default 
        #when the query is first run without the output_file argument, and the output file 
        #path is saved in the same directory as the pid file.
        $out_file = $output_file // "$tmpdir/output.html.info.txt";
        
        #Note: the temporary sas file is named "tmp_sas_codes_{timestamp}.sas" when the sas codes
        # contain special characters and need to be saved into a temporary sas file to run the query, 
        #and the temporary sas file is saved in the same directory as the pid file.
        $tmp_sas_file = $args->{tmp_sas_file} // "$tmpdir/tmp_sas_file.sas"; # This is for storing the temporary sas file path when the sas codes contain special characters;
        
        unless (-f $pid_file) {
            my $pid = fork();
            return { content => [{ type => "text", text => "ERROR: Could not fork." }] } unless defined $pid;

            if ($pid == 0) {
                local $ENV{PIPELINE_SAS_ODA_ACCOUNT} = $sas_oda_account if defined $sas_oda_account && length $sas_oda_account;
                local $ENV{PIPELINE_SAS_ODA_PASSWORD} = $sas_oda_password if defined $sas_oda_password && length $sas_oda_password;
                local $ENV{PIPELINE_FORCE_SAS_ODA_AUTH_PROMPT} = 1 if defined $prompt_sas_oda_auth && $prompt_sas_oda_auth =~ /^(?:1|true|yes|y)$/i;
                my $cmd = "perl -S run_sas_codes_or_script_in_ODA.pl --output-prefix $tmpdir ";
                $cmd .= "--persistent --session-id " . shell_quote_single($session_id_arg) . " " if $use_persistent;
                $cmd .= $sas_codes;
                print STDERR "Executing command in background: $cmd\n";
                exec "$cmd 2>&1";
                exit(1);
            }

            open(my $pfh, '>', $pid_file);
            print $pfh $pid;
            close $pfh;
            # Return initial status with PID and instructions for checking later, avoiding the endless loop of checking by the AI agent, 
            # which may cause excessive checking and potential issues with the server;
            # If the sas_codes_or_file contains single/double quotes, the tool will save the sas codes into a temporary sas file;
            if (-f $tmp_sas_file){
              return {
                content => [{
                    type => "text",
                    text => "QUERYING: SAS ODA for $sas_codes\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"sas_codes\": \"$tmp_sas_file\", \"pid\": $pid}\n" . 
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"
                }]
            };
            }else{
              return {
                content => [{
                    type => "text",
                    text => "QUERYING: SAS ODA for $sas_codes_or_file\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"sas_codes_or_file\": \"$sas_codes_or_file\", \"pid\": $pid}\n" .
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"}]
        };
       }
    }
  }
);

$server->tool(
    name        => 'auto_prepare_and_run_diff_gwas',
    description => 'Run the generalized differential-GWAS automation pipeline from a comparison spec JSON file or from an auto-detected GWAS directory. ' .
    'This wraps auto_prepare_and_run_diff_gwas.pl as a background MCP job so users can prepare configs and run the ' .
    'GWAS comparison workflow, including SAS ODA plots, by providing either a spec file or a gwas_dir for auto-detection. ' .
    'The main entrypoint stays at the project root, while its helper scripts are organized under DiffGWASDeps/ for easier maintenance. ' .
    'It now also supports merged-wide GWAS tables through source_mode=merged_gwas_table, where one shared-locus table contains cohort-level BETA/SE/P blocks such as BETA_DS_ALL, SE_DS_ALL, P_DS_ALL and BETA_MP2PRT, SE_MP2PRT, P_MP2PRT, plus optional extra association tracks such as meta-analysis P/Z columns. ' .
    'That merged-wide mode can preview or write an inferred spec JSON from gwas_dir, normalize the merged table into the plotting-wide schema used by the rest of the repository, and then drive genome-wide Manhattan, local Manhattan, local GTF, and forest plots from the same automation entrypoint. ' .
    'The same workflow now supports full multi-GWAS comparison views, single selected GWAS association views, and inquiry-SNP local Manhattan / local GTF views through the shared display_gwas and target_snps controls. ' .
    'For custom genome-wide Manhattan subsets, the plotting path expects chromosome labels to be normalized into sortable numeric coordinates before SAS reads them; leaving literal X/Y values in a numeric CHR import can create a false pre-chr1 block in the final figure. ' .
    'The local GTF path now also understands reference genome build selection through an explicit reference_build override or automatic header/path-token detection, and then maps hg19, hg38, or T2T/hs1 inputs onto the matching built-in GTF profile. ' .
    'It supports both differential top-hit mode and common-association top-hit mode, where loci can be selected from strong single-GWAS associations that also have nominal association in another GWAS. ' .
    'It can also forward optional SNP:GENE label overrides and emit runnable local desktop SAS scripts instead of submitting everything to SAS ODA when users prefer workstation SAS. ' .
    'It can also forward local-Manhattan layout controls so users can manually tune crowded SNP-gene labels without editing SAS templates by hand. ' .
    'For local GTF plots, it can forward a dedicated local GTF window override that now controls both the extracted GTF subset and the displayed local gene-track window. ' .
    'The batched and single-SNP SAS local-GTF wrappers now also share the same more readable baseline figure defaults, starting from GTF_PCT4NEG_Y=1.4 and a practical 1000-pixel local-GTF design height before any per-locus auto-tuning. ' .
    'Large-window local GTF reruns now use a gzip-compressed pre-extracted GTF subset for SAS ODA upload, which is much more reliable than constructing the full region inside SAS WORK. ' .
    'The local GTF path now also prefers the repository copy of SNP_Local_Manhattan_With_GTF.sas and forces the final displayed x-axis back to the observed association-signal span, which helps avoid chr plots that incorrectly start at 0 after aggressive gene-track expansion. ' .
    'If a long local GTF submit finishes remotely but the expected HTML download path is flaky, the wrapper can now recover from the helper-saved sas_res_*.html artifact and the PNG path reported in the SAS log. When that PNG exists, the opened result now prefers a small figure-first final HTML wrapper and keeps the raw SAS HTML as a sidecar. ' .
    'On the first SAS-backed run, the vendored SAS ODA helper now prompts for missing credentials or accepts optional account/password fields, validates them with proc setinit;run;, and saves the successful login for later reuse. ' .
    'The same automation path has now also been validated from an isolated Ubuntu Docker runtime; for Linux local-GTF runs, keeping cache/gtf/gencode.v49.annotation.gtf.gz available locally avoids treating a transient EBI download failure as a pipeline failure. ' .
    'In Ubuntu Docker, the first image build can take several minutes because the container provisions system packages plus repo-local Perl/Python runtimes, and long SAS ODA local-GTF reruns can also spend extra time in remote housekeeping after the scientific plot itself has already finished. ' .
    'Users who prefer a containerized Linux deployment can now start from the bundled Dockerfile or install/singularity/MultiGWAS-Explorer_pipeline.def instead of reconstructing the runtime manually. ' .
    'It now also supports Nextflow-like targeted reruns: list available step names, run only one named stage, or rerun a step range without recomputing the whole workflow. ' .
    'The repo now prefers vendored helper scripts plus repo-local Perl and Python runtimes when present, including platform-specific local Perl trees such as local/perl5-cygwin and local/perl5-linux so Windows portable Cygwin does not accidentally load Linux GD or zlib modules. ' .
    'When one automation call requests multiple plot stages that consume the same wide gz subset, the pipeline now uploads that data file to SAS ODA once, reuses it across the requested plot runners, and deletes it once at the end if cleanup is enabled. ' .
    'The same tool can be called again with a PID to check status and retrieve the final log/output summary.',
    input_schema => {
        type => 'object',
        properties => {
            spec_file => {
                type => 'string',
                description => 'Optional comparison spec JSON file path. The spec describes source_mode, input/output locations, groups, pairs, and optional plotting settings.'
            },
            gwas_dir => {
                type => 'string',
                description => 'Optional GWAS directory to scan and auto-detect. This is the preferred entrypoint for merged-wide GWAS tables such as the AOA merged DS_ALL + MP2PRT + meta table.'
            },
            spec_out => {
                type => 'string',
                description => 'Optional output path for an auto-generated spec JSON when using gwas_dir, for example: configs/auto_aoa_merged.spec.json'
            },
            raw_column_alias_config => {
                type => 'string',
                description => 'Optional JSON file that adds raw header aliases for new GWAS formats during auto-detection.'
            },
            generate_spec_only => {
                type => 'string',
                description => 'Optional truthy flag to stop after generating or previewing the inferred spec. This is especially useful for merged-wide GWAS tables when users want to inspect the inferred source_mode, cohort blocks, and output paths before running plots.'
            },
            preview_spec => {
                type => 'string',
                description => 'Optional truthy flag to print the inferred spec JSON to the log or response before writing it. Recommended first step for merged-wide AOA-style inputs.'
            },
            project_tag => {
                type => 'string',
                description => 'Optional override for the inferred project_tag during gwas_dir auto-detection.'
            },
            artifact_stem => {
                type => 'string',
                description => 'Optional override for the inferred artifact_stem during gwas_dir auto-detection.'
            },
            mode => {
                type => 'string',
                description => 'Optional mode for auto_prepare_and_run_diff_gwas.pl. Typical values: full or configs. Default: full.'
            },
            plots => {
                type => 'string',
                description => 'Optional comma-separated plot list, for example: manhattan,local_manhattan,local_gtf,forest'
            },
            skip_plots => {
                type => 'string',
                description => 'Optional truthy flag to skip SAS ODA plotting. Accepts values like 1, true, yes.'
            },
            force => {
                type => 'string',
                description => 'Optional truthy flag to rerun steps even when expected outputs already exist.'
            },
            list_steps => {
                type => 'string',
                description => 'Optional truthy flag to print the available step names for the current spec and exit. Good for helping the LLM decide which exact stage to rerun.'
            },
            step => {
                type => 'string',
                description => 'Optional step name or comma-separated step names to run selectively, for example: extract_wide_subset or plot_local_gtf. This is the most direct way to rerun one specific stage.'
            },
            from_step => {
                type => 'string',
                description => 'Optional first step in a rerun range, for example: extract_wide_subset'
            },
            to_step => {
                type => 'string',
                description => 'Optional last step in a rerun range, for example: plot_local_gtf'
            },
            merge_raw => {
                type => 'string',
                description => 'Optional truthy convenience flag for --merge-raw'
            },
            sort_long => {
                type => 'string',
                description => 'Optional truthy convenience flag for --sort-long'
            },
            diff_pairs => {
                type => 'string',
                description => 'Optional truthy convenience flag for --diff-pairs'
            },
            standardize_diff => {
                type => 'string',
                description => 'Optional truthy convenience flag for --standardize-diff'
            },
            extract_wide_subset => {
                type => 'string',
                description => 'Optional truthy convenience flag for --extract-wide-subset'
            },
            plot_manhattan => {
                type => 'string',
                description => 'Optional truthy convenience flag for --plot-manhattan'
            },
            plot_local_manhattan => {
                type => 'string',
                description => 'Optional truthy convenience flag for --plot-local-manhattan'
            },
            plot_local_gtf => {
                type => 'string',
                description => 'Optional truthy convenience flag for --plot-local-gtf'
            },
            plot_forest => {
                type => 'string',
                description => 'Optional truthy convenience flag for --plot-forest'
            },
            local_max_hits_per_fig => {
                type => 'integer',
                description => 'Optional requested upper bound for local top-hit columns per panel. The current pipeline allows up to 30 columns per local Manhattan figure.'
            },
            reference_build => {
                type => 'string',
                description => 'Optional reference genome build override for build-aware local GTF annotation, for example: hg19, hg38, or t2t. If omitted, the pipeline tries explicit spec fields first, then header/path tokens, and otherwise falls back to hg38.'
            },
            local_gtf_window_bp => {
                type => 'string',
                description => 'Optional override for the half-window size used by local GTF plots, for example: 1e8 or 2e7. This can differ from the local Manhattan window, and now controls both the extracted local GTF subset and the displayed local GTF plot range. Large-window reruns use a gzip-compressed local GTF subset upload to SAS ODA, while the current SAS wrappers keep a shared larger bottom gene-track baseline for readability.'
            },
            local_manhattan_angle4xaxis_label => {
                type => 'string',
                description => 'Optional override for local Manhattan SNP-gene x-axis label rotation, for example: 60.'
            },
            local_manhattan_xgrp_y_pos => {
                type => 'string',
                description => 'Optional override for the vertical position of local Manhattan SNP-gene labels, for example: -2.5.'
            },
            local_manhattan_yoffset_top => {
                type => 'string',
                description => 'Optional override for the upper local Manhattan label offset.'
            },
            local_manhattan_yoffset_bottom => {
                type => 'string',
                description => 'Optional override for the lower local Manhattan label offset.'
            },
            local_manhattan_fontsize => {
                type => 'string',
                description => 'Optional override for local Manhattan SNP-gene label font size.'
            },
            local_manhattan_y_axis_label_size => {
                type => 'string',
                description => 'Optional override for local Manhattan y-axis title size.'
            },
            local_manhattan_y_axis_value_size => {
                type => 'string',
                description => 'Optional override for local Manhattan y-axis tick-label size.'
            },
            get_common_associations => {
                type => 'string',
                description => 'Optional flag or threshold for common-association local-hit mode. Use values like true, yes, or 5e-8. This selects loci with strong single-GWAS association plus nominal association in another GWAS with the same effect direction.'
            },
            common_association_top_hit_threshold => {
                type => 'string',
                description => 'Optional explicit starting threshold for common-association local-hit mode, for example: 5e-8'
            },
            display_gwas => {
                type => 'string',
                description => 'Optional comma-separated GWAS track selection shared by the SAS ODA and gunplot pipelines. Use pair prefixes such as ALL,EUR,ASN for differential tracks and GWAS labels such as ALL_FEMALE or EUR_MALE for single-GWAS tracks.'
            },
            target_snps => {
                type => 'string',
                description => 'Optional comma-separated inquiry SNP list to drive local Manhattan and local GTF plots directly, for example: rs17425819,rs185665940'
            },
            target_snp_genes => {
                type => 'string',
                description => 'Optional comma-separated SNP:GENE overrides, for example: rs17425819:JAK2,rs185665940:FANCL'
            },
            sas_oda_account => {
                type => 'string',
                description => 'Optional SAS ODA account/email for noninteractive first-run login bootstrap.'
            },
            sas_oda_password => {
                type => 'string',
                description => 'Optional SAS ODA password for noninteractive first-run login bootstrap.'
            },
            prompt_sas_oda_auth => {
                type => 'string',
                description => 'Optional truthy flag to force a SAS ODA credential refresh even when a saved authinfo entry already exists.'
            },
            emit_local_sas_scripts => {
                type => 'string',
                description => 'Optional truthy flag to emit local desktop-SAS runnable plot scripts alongside the SAS ODA workflow.'
            },
            local_sas_only => {
                type => 'string',
                description => 'Optional truthy flag to emit local desktop-SAS runnable plot scripts and stop before SAS ODA upload/submit work.'
            },
            exclude_non_protein_coding_genes_in_local_gtf => {
                type => 'string',
                description => 'Legacy truthy convenience flag for --exclude-non-protein-coding-genes-in-local-gtf. Local GTF plots are protein-coding-only by default; set include_non_protein_coding_genes_in_local_gtf=1 in the spec JSON when non-coding genes should be included.'
            },
            cleanup_shared_plot_data => {
                type => 'string',
                description => 'Optional truthy convenience flag for --cleanup-shared-plot-data'
            },
            output_file => {
                type => 'string',
                description => 'Optional output log file path. By default a tmp*/output.html.info.txt file is created.'
            },
            pid => {
                type => 'integer',
                description => 'Optional PID to check status of a previous auto_prepare_and_run_diff_gwas query.'
            }
        },
        required => ['spec_file']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $spec_file = $args->{spec_file} // '';
        my $gwas_dir = $args->{gwas_dir} // '';
        my $spec_out = $args->{spec_out} // '';
        my $raw_column_alias_config = $args->{raw_column_alias_config} // '';
        my $generate_spec_only = $args->{generate_spec_only} // '';
        my $preview_spec = $args->{preview_spec} // '';
        my $project_tag = $args->{project_tag} // '';
        my $artifact_stem = $args->{artifact_stem} // '';
        my $mode = $args->{mode} // 'full';
        my $plots = $args->{plots} // '';
        my $skip_plots = $args->{skip_plots} // '';
        my $force = $args->{force} // '';
        my $list_steps = $args->{list_steps} // '';
        my $step = $args->{step} // '';
        my $from_step = $args->{from_step} // '';
        my $to_step = $args->{to_step} // '';
        my $local_max_hits_per_fig = $args->{local_max_hits_per_fig};
        my $reference_build = $args->{reference_build} // '';
        my $local_gtf_window_bp = $args->{local_gtf_window_bp} // '';
        my $local_manhattan_angle4xaxis_label = $args->{local_manhattan_angle4xaxis_label} // '';
        my $local_manhattan_xgrp_y_pos = $args->{local_manhattan_xgrp_y_pos} // '';
        my $local_manhattan_yoffset_top = $args->{local_manhattan_yoffset_top} // '';
        my $local_manhattan_yoffset_bottom = $args->{local_manhattan_yoffset_bottom} // '';
        my $local_manhattan_fontsize = $args->{local_manhattan_fontsize} // '';
        my $local_manhattan_y_axis_label_size = $args->{local_manhattan_y_axis_label_size} // '';
        my $local_manhattan_y_axis_value_size = $args->{local_manhattan_y_axis_value_size} // '';
        my $get_common_associations = $args->{get_common_associations} // '';
        my $common_association_top_hit_threshold = $args->{common_association_top_hit_threshold} // '';
        my $display_gwas = $args->{display_gwas} // '';
        my $target_snps = $args->{target_snps} // '';
        my $target_snp_genes = $args->{target_snp_genes} // '';
        my $sas_oda_account = $args->{sas_oda_account} // '';
        my $sas_oda_password = $args->{sas_oda_password} // '';
        my $prompt_sas_oda_auth = $args->{prompt_sas_oda_auth} // '';
        my $emit_local_sas_scripts = $args->{emit_local_sas_scripts} // '';
        my $local_sas_only = $args->{local_sas_only} // '';
        my $exclude_non_protein_coding_genes_in_local_gtf = $args->{exclude_non_protein_coding_genes_in_local_gtf} // '';
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        my $out_file;
        my $pid_file;

        if (defined $pid_arg) {
            $pid_file = "tmp*/auto_prepare_and_run_diff_gwas.pid";
            my @pid_files = glob($pid_file);

            foreach my $pf (@pid_files) {
                next unless -f $pf;
                open(my $pfh, '<', $pf);
                my $stored_line = <$pfh>;
                close $pfh;
                chomp $stored_line;

                my ($stored_pid, $stored_out_file, $stored_err_file) = split(/\t/, $stored_line, 3);
                $stored_out_file //= '';
                $stored_err_file //= '';

                next unless $stored_pid == $pid_arg;

                if (defined $stored_out_file && length $stored_out_file) {
                    $out_file = $stored_out_file;
                } else {
                    $out_file = $pf;
                    $out_file =~ s/\.pid$/.txt/;
                }

                if (pid_is_running($pid_arg)) {
                    return {
                        content => [{
                            type => "text",
                            text => "STATUS: RUNNING (PID $pid_arg)\n" .
                                    "Output file: $out_file\nAsk the AI agent to check status again in a moment."
                        }]
                    };
                }

                my $content = "";
                if (-f $out_file && open(my $fh, '<', $out_file)) {
                    local $/;
                    $content = <$fh>;
                    close $fh;
                }
                if (defined $stored_err_file && length $stored_err_file && -f $stored_err_file) {
                    if (open(my $efh, '<', $stored_err_file)) {
                        local $/;
                        my $err_content = <$efh>;
                        close $efh;
                        if (defined $err_content && length $err_content) {
                            if (length $content) {
                                $content .= "\n\n[stderr]\n" . $err_content;
                            }
                            else {
                                $content = $err_content;
                            }
                        }
                    }
                }

                unlink $pf;
                cleanup_generated_tmpdir_if_empty($pf);

                return {
                    content => [{
                        type => "text",
                        text => "STATUS: COMPLETE (PID $pid_arg)\n\nAutomation log saved to: $out_file\n\n" . $content
                    }]
                };
            }

            return {
                content => [{
                    type => "text",
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        die "Either spec_file or gwas_dir is required\n" unless length($spec_file) || length($gwas_dir);

        my $tmpdir = "tmp" . time();
        unless (-d $tmpdir) {
            mkdir $tmpdir or die "Failed to create tmp directory: $!";
        }

        $pid_file = "$tmpdir/auto_prepare_and_run_diff_gwas.pid";
        $out_file = $output_file // "$tmpdir/output.html.info.txt";

        my $out_dir = dirname($out_file);
        if (defined $out_dir && length $out_dir && $out_dir ne '.' && !-d $out_dir) {
            make_path($out_dir) or die "Failed to create output directory $out_dir: $!";
        }

        unless (-f $pid_file) {
            my $automation_script = File::Spec->catfile($Bin, 'auto_prepare_and_run_diff_gwas.pl');
            my @cmd = (
                'perl',
                $automation_script,
            );
            push @cmd, ('--spec', $spec_file) if defined $spec_file && length $spec_file;
            push @cmd, ('--gwas-dir', $gwas_dir) if defined $gwas_dir && length $gwas_dir;
            push @cmd, ('--spec-out', $spec_out) if defined $spec_out && length $spec_out;
            push @cmd, ('--raw-column-alias-config', $raw_column_alias_config)
              if defined $raw_column_alias_config && length $raw_column_alias_config;
            push @cmd, '--generate-spec-only'
              if defined $generate_spec_only && $generate_spec_only =~ /^(?:1|true|yes|y)$/i;
            push @cmd, '--preview-spec'
              if defined $preview_spec && $preview_spec =~ /^(?:1|true|yes|y)$/i;
            push @cmd, ('--project-tag', $project_tag) if defined $project_tag && length $project_tag;
            push @cmd, ('--artifact-stem', $artifact_stem) if defined $artifact_stem && length $artifact_stem;
            push @cmd, ('--mode', $mode) if defined $mode && length $mode;
            push @cmd, ('--plots', $plots) if defined $plots && length $plots;
            push @cmd, '--skip-plots' if defined $skip_plots && $skip_plots =~ /^(?:1|true|yes|y)$/i;
            push @cmd, '--force' if defined $force && $force =~ /^(?:1|true|yes|y)$/i;
            push @cmd, '--list-steps' if defined $list_steps && $list_steps =~ /^(?:1|true|yes|y)$/i;
            push @cmd, ('--step', $step) if defined $step && length $step;
            push @cmd, ('--from-step', $from_step) if defined $from_step && length $from_step;
            push @cmd, ('--to-step', $to_step) if defined $to_step && length $to_step;
            push @cmd, ('--local-max-hits-per-fig', $local_max_hits_per_fig) if defined $local_max_hits_per_fig && $local_max_hits_per_fig =~ /^\d+$/;
            push @cmd, ('--reference-build', $reference_build)
              if defined $reference_build && length $reference_build;
            push @cmd, ('--local-gtf-window-bp', $local_gtf_window_bp)
              if defined $local_gtf_window_bp && length $local_gtf_window_bp;
            push @cmd, ('--local-manhattan-angle4xaxis-label', $local_manhattan_angle4xaxis_label)
              if defined $local_manhattan_angle4xaxis_label && length $local_manhattan_angle4xaxis_label;
            push @cmd, ('--local-manhattan-xgrp-y-pos', $local_manhattan_xgrp_y_pos)
              if defined $local_manhattan_xgrp_y_pos && length $local_manhattan_xgrp_y_pos;
            push @cmd, ('--local-manhattan-yoffset-top', $local_manhattan_yoffset_top)
              if defined $local_manhattan_yoffset_top && length $local_manhattan_yoffset_top;
            push @cmd, ('--local-manhattan-yoffset-bottom', $local_manhattan_yoffset_bottom)
              if defined $local_manhattan_yoffset_bottom && length $local_manhattan_yoffset_bottom;
            push @cmd, ('--local-manhattan-fontsize', $local_manhattan_fontsize)
              if defined $local_manhattan_fontsize && length $local_manhattan_fontsize;
            push @cmd, ('--local-manhattan-y-axis-label-size', $local_manhattan_y_axis_label_size)
              if defined $local_manhattan_y_axis_label_size && length $local_manhattan_y_axis_label_size;
            push @cmd, ('--local-manhattan-y-axis-value-size', $local_manhattan_y_axis_value_size)
              if defined $local_manhattan_y_axis_value_size && length $local_manhattan_y_axis_value_size;
            if (defined $get_common_associations && length $get_common_associations) {
                if ($get_common_associations =~ /^(?:1|true|yes|y)$/i) {
                    push @cmd, '--get-common-associations';
                } else {
                    push @cmd, "--get-common-associations=$get_common_associations";
                }
            }
            push @cmd, ('--common-association-top-hit-threshold', $common_association_top_hit_threshold)
              if defined $common_association_top_hit_threshold && length $common_association_top_hit_threshold;
            push @cmd, ('--display-gwas', $display_gwas)
              if defined $display_gwas && length $display_gwas;
            push @cmd, ('--target-snps', $target_snps)
              if defined $target_snps && length $target_snps;
            push @cmd, ('--target-snp-genes', $target_snp_genes)
              if defined $target_snp_genes && length $target_snp_genes;
            push @cmd, '--emit-local-sas-scripts'
              if defined $emit_local_sas_scripts && $emit_local_sas_scripts =~ /^(?:1|true|yes|y)$/i;
            push @cmd, '--local-sas-only'
              if defined $local_sas_only && $local_sas_only =~ /^(?:1|true|yes|y)$/i;
            push @cmd, '--exclude-non-protein-coding-genes-in-local-gtf'
              if defined $exclude_non_protein_coding_genes_in_local_gtf
                 && $exclude_non_protein_coding_genes_in_local_gtf =~ /^(?:1|true|yes|y)$/i;
            for my $step_flag (
                [merge_raw => '--merge-raw'],
                [sort_long => '--sort-long'],
                [diff_pairs => '--diff-pairs'],
                [standardize_diff => '--standardize-diff'],
                [extract_wide_subset => '--extract-wide-subset'],
                [plot_manhattan => '--plot-manhattan'],
                [plot_local_manhattan => '--plot-local-manhattan'],
                [plot_local_gtf => '--plot-local-gtf'],
                [plot_forest => '--plot-forest'],
                [cleanup_shared_plot_data => '--cleanup-shared-plot-data'],
            ) {
                my ($arg_key, $cli_flag) = @{$step_flag};
                my $val = $args->{$arg_key} // '';
                push @cmd, $cli_flag if defined $val && $val =~ /^(?:1|true|yes|y)$/i;
            }

            my $cmd_preview = join(' ', map { shell_quote_single($_) } @cmd);
            my $pid;

            if ($^O =~ /^(?:cygwin|MSWin32)$/i) {
                my $perl_exe_win = cygpath_to_windows($^X);
                my $cwd_win = cygpath_to_windows(getcwd());
                my $out_file_abs = File::Spec->rel2abs($out_file);
                my $err_file = $out_file . '.stderr.log';
                my $err_file_abs = File::Spec->rel2abs($err_file);
                my $out_file_win = cygpath_to_windows($out_file_abs);
                my $err_file_win = cygpath_to_windows($err_file_abs);
                my $launch_pid_file = File::Spec->catfile($tmpdir, 'auto_prepare_launch.pid');
                my $launch_pid_file_win = cygpath_to_windows(File::Spec->rel2abs($launch_pid_file));
                my $arg_list = '@(' . join(', ', map { powershell_quote_single($_) } @cmd[1 .. $#cmd]) . ')';
                my @ps_env_prefix;
                push @ps_env_prefix, '$env:PIPELINE_SAS_ODA_ACCOUNT=' . powershell_quote_single($sas_oda_account) . ';'
                  if defined $sas_oda_account && length $sas_oda_account;
                push @ps_env_prefix, '$env:PIPELINE_SAS_ODA_PASSWORD=' . powershell_quote_single($sas_oda_password) . ';'
                  if defined $sas_oda_password && length $sas_oda_password;
                push @ps_env_prefix, '$env:PIPELINE_FORCE_SAS_ODA_AUTH_PROMPT=1;'
                  if defined $prompt_sas_oda_auth && $prompt_sas_oda_auth =~ /^(?:1|true|yes|y)$/i;
                my $ps_cmd = join ' ',
                    @ps_env_prefix,
                    '$p = Start-Process',
                    '-FilePath', powershell_quote_single($perl_exe_win),
                    '-ArgumentList', $arg_list,
                    '-WorkingDirectory', powershell_quote_single($cwd_win),
                    '-RedirectStandardOutput', powershell_quote_single($out_file_win),
                    '-RedirectStandardError', powershell_quote_single($err_file_win),
                    '-WindowStyle Hidden -PassThru;',
                    'Set-Content -Path', powershell_quote_single($launch_pid_file_win), '-Value $p.Id;';
                print STDERR "Executing command in background via Start-Process: $cmd_preview\n";
                my $ps_status = system('powershell', '-NoProfile', '-NonInteractive', '-Command', $ps_cmd);
                die "Failed to invoke PowerShell Start-Process launcher\n" if $ps_status != 0;
                open(my $pidfh, '<', $launch_pid_file)
                  or die "Failed to read launched background PID from $launch_pid_file: $!\n";
                my $pid_out = <$pidfh>;
                close $pidfh;
                unlink $launch_pid_file;
                chomp $pid_out;
                ($pid) = ($pid_out =~ /(\d+)/);
                die "Failed to start background automation process via PowerShell Start-Process\n" unless $pid;
                open(my $pfh, '>', $pid_file);
                print $pfh "$pid\t$out_file\t$err_file\n";
                close $pfh;
            }
            else {
                my $child_pid = fork();
                return { content => [{ type => "text", text => "ERROR: Could not fork." }] } unless defined $child_pid;

                if ($child_pid == 0) {
                    local $ENV{PIPELINE_SAS_ODA_ACCOUNT} = $sas_oda_account if defined $sas_oda_account && length $sas_oda_account;
                    local $ENV{PIPELINE_SAS_ODA_PASSWORD} = $sas_oda_password if defined $sas_oda_password && length $sas_oda_password;
                    local $ENV{PIPELINE_FORCE_SAS_ODA_AUTH_PROMPT} = 1 if defined $prompt_sas_oda_auth && $prompt_sas_oda_auth =~ /^(?:1|true|yes|y)$/i;
                    open(my $child_out, '>', $out_file) or die "Failed to open output file $out_file: $!";
                    open(STDIN, '<', File::Spec->devnull()) or die "Failed to redirect STDIN to devnull: $!";
                    open(STDOUT, '>&', $child_out) or die "Failed to redirect STDOUT to $out_file: $!";
                    open(STDERR, '>&', \*STDOUT) or die "Failed to redirect STDERR to STDOUT for $out_file: $!";
                    select STDOUT; $| = 1;
                    select STDERR; $| = 1;
                    print STDERR "Executing command in background: $cmd_preview\n";
                    exec @cmd;
                    exit(1);
                }
                $pid = $child_pid;
                open(my $pfh, '>', $pid_file);
                print $pfh "$pid\t$out_file\n";
                close $pfh;
            }

            return {
                content => [{
                        type => "text",
                        text => "QUERYING: auto_prepare_and_run_diff_gwas for " . (length($spec_file) ? $spec_file : $gwas_dir) . "\nPID: $pid\nOutput file: $out_file\n" .
                            "Ask the AI agent to check status with: {" .
                            (length($spec_file) ? "\"spec_file\": \"$spec_file\", " : '') .
                            (length($gwas_dir) ? "\"gwas_dir\": \"$gwas_dir\", " : '') .
                            "\"pid\": $pid}\n" .
                            "Hints: use preview_spec=true and generate_spec_only=true with gwas_dir for merged-wide GWAS inspection; use spec_out=\"configs/auto_aoa_merged.spec.json\" to keep the inferred spec; use list_steps=true to inspect exact stage names; use step=\"plot_local_gtf\" or from_step=\"extract_wide_subset\" for targeted reruns; use display_gwas=\"ALL_FEMALE\" or display_gwas=\"EUR,EUR_FEMALE,EUR_MALE\" for custom displayed GWAS tracks; use target_snps=\"rs123\" for inquiry-SNP local plots; use get_common_associations=\"true\" or get_common_associations=\"5e-8\" for shared-association top-hit mode with concordant direction.\n" .
                            "The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking.\n"
                }]
            };
        }
    }
);


# # -------------------------
# # Run the gunplot wrapper (PDL + gnuplot) as an asynchronous MCP job
# # -------------------------
$server->tool(
    name        => 'run_gunplot_wrapper',
    description => 'Run the auto_prepare_and_run_diff_gwas_with_gunplot.pl wrapper to produce genomewide Manhattan, local Manhattan, and local GTF plots with the alternative gunplot renderer. This backend now mirrors the SAS ODA path for displayed-GWAS selection, single-GWAS rendering mode, inquiry SNP local plots, optional SNP:GENE label overrides, and merged-wide GWAS compatibility once a spec has been generated. For merged-wide study tables such as the AOA DS_ALL + MP2PRT + meta input, first generate the inferred spec through auto_prepare_and_run_diff_gwas.pl with gwas_dir/preview_spec/spec_out, then pass that spec here. The genomewide Manhattan renderer now uses the same repeated chromosome palette family and top-of-panel GWAS labels as the SAS ODA multi-track figure style, while still retaining small renderer-specific differences such as gnuplot rasterization. It also shares the build-aware local GTF logic, so explicit reference_build overrides or detected hg19/hg38/T2T tokens select the matching built-in GTF profile automatically. The wrapper now prefers a real gnuplot command from the active PATH, keeps portable Cygwin installs self-contained instead of borrowing a host-specific Windows gnuplot.exe, validates cached genome-wide wide subsets against their manifest before reuse, has been exercised successfully from an isolated Ubuntu Docker runtime as well as portable Cygwin, and can also be packaged through the bundled Dockerfile or Singularity/Apptainer definition. In Ubuntu Docker, the first image build can take several minutes, and the genomewide Manhattan stage is usually the dominant runtime cost; for quick Linux-container validation, prefer one inquiry SNP with local_manhattan and local_gtf before adding manhattan.',
    input_schema => {
        type => 'object',
        properties => {
            spec_file => { type => 'string', description => 'Required comparison spec JSON file path. For merged-wide GWAS tables, generate this spec first through auto_prepare_and_run_diff_gwas.pl --gwas-dir ... --spec-out ... .' },
            plots => { type => 'string', description => 'Optional comma-separated plot list, for example: manhattan,local_manhattan,local_gtf' },
            step => { type => 'string', description => 'Optional plot step name or comma-separated step names, for example: plot_manhattan or plot_local_gtf' },
            force => { type => 'string', description => 'Optional truthy flag to refresh cached preprocessing and rerender outputs.' },
            remove_x_chr => { type => 'string', description => 'Optional truthy flag to remove chrX from final gunplot figures. Default: true.' },
            display_gwas => { type => 'string', description => 'Optional comma-separated displayed GWAS tracks. Use pair prefixes such as ALL,EUR,ASN for differential tracks and GWAS labels such as ALL_FEMALE or EUR_MALE for single-GWAS tracks.' },
            target_snps => { type => 'string', description => 'Optional comma-separated inquiry SNP list for local Manhattan / local GTF plots. This is also the preferred fast-validation path for Ubuntu Docker runs.' },
            target_snp_genes => { type => 'string', description => 'Optional comma-separated SNP:GENE overrides, for example: rs17425819:JAK2,rs185665940:FANCL' },
            get_common_associations => { type => 'string', description => 'Optional flag or threshold for common-association local-hit mode, for example: true or 5e-8.' },
            common_association_top_hit_threshold => { type => 'string', description => 'Optional explicit starting threshold for common-association local-hit mode.' },
            reference_build => { type => 'string', description => 'Optional reference genome build override for build-aware local GTF annotation, for example: hg19, hg38, or t2t. If omitted, the pipeline tries explicit spec fields first, then header/path tokens, and otherwise falls back to hg38.' },
            local_gtf_window_bp => { type => 'string', description => 'Optional override for local GTF half-window size.' },
            local_max_hits_per_fig => { type => 'integer', description => 'Optional override for local Manhattan batching.' },
            local_manhattan_columns => { type => 'integer', description => 'Optional override for combined local-Manhattan figure columns.' },
            local_manhattan_annotation => { type => 'string', description => 'Optional under-column annotation mode for combined local Manhattan: labels, gtf, auto, none.' },
            output_file => { type => 'string', description => 'Optional path to write wrapper stdout/stderr.' },
            pid       => { type => 'integer', description => 'Optional PID to check status of a previous gunplot wrapper query.' }
        },
        required => []
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $spec = $args->{spec_file} // '';
        my $plots = $args->{plots} // '';
        my $step = $args->{step} // '';
        my $force = $args->{force} // '';
        my $remove_x_chr = $args->{remove_x_chr};
        my $display_gwas = $args->{display_gwas} // '';
        my $target_snps = $args->{target_snps} // '';
        my $target_snp_genes = $args->{target_snp_genes} // '';
        my $get_common_associations = $args->{get_common_associations} // '';
        my $common_association_top_hit_threshold = $args->{common_association_top_hit_threshold} // '';
        my $reference_build = $args->{reference_build} // '';
        my $local_gtf_window_bp = $args->{local_gtf_window_bp} // '';
        my $local_max_hits_per_fig = $args->{local_max_hits_per_fig};
        my $local_manhattan_columns = $args->{local_manhattan_columns};
        my $local_manhattan_annotation = $args->{local_manhattan_annotation} // '';
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        my $out_file;
        my $pid_file;

        if (defined $pid_arg) {
            $pid_file = "tmp*/gunplot_wrapper.pid";
            my @pid_files = glob($pid_file);
            foreach my $pf (@pid_files) {
                next unless -f $pf;
                open(my $pfh, '<', $pf);
                my $stored_line = <$pfh>;
                close $pfh;
                chomp $stored_line;
                my ($stored_pid, $stored_out_file) = split(/\t/, $stored_line, 2);
                next unless $stored_pid == $pid_arg;
                $out_file = $stored_out_file;
                if (pid_is_running($pid_arg)) {
                    return {
                        content => [{
                            type => 'text',
                            text => "STATUS: RUNNING (PID $pid_arg)\nOutput file: $out_file\nAsk the AI agent to check status again in a moment."
                        }]
                    };
                }
                my $content = '';
                if (-f $out_file && open(my $fh, '<', $out_file)) {
                    local $/;
                    $content = <$fh>;
                    close $fh;
                }
                unlink $pf;
                cleanup_generated_tmpdir_if_empty($pf);
                return {
                    content => [{
                        type => 'text',
                        text => "STATUS: COMPLETE (PID $pid_arg)\n\nGunplot log saved to: $out_file\n\n" . $content
                    }]
                };
            }
            return {
                content => [{
                    type => 'text',
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        die "spec_file is required\n" unless length($spec);

        my $tmpdir = "tmp" . time();
        mkdir $tmpdir or die "Failed to create tmp dir: $!" unless -d $tmpdir;
        $pid_file = "$tmpdir/gunplot_wrapper.pid";
        $out_file = $output_file // "$tmpdir/output.gunplot.txt";

        my $script = File::Spec->catfile($Bin, 'auto_prepare_and_run_diff_gwas_with_gunplot.pl');
        my @cmd = ('perl', $script, '--spec', $spec);
        push @cmd, ('--plots', $plots) if length $plots;
        if (length $step) {
            for my $step_name (grep { length } split /\s*,\s*/, $step) {
                push @cmd, ('--step', $step_name);
            }
        }
        push @cmd, '--force' if defined $force && $force =~ /^(?:1|true|yes|y)$/i;
        if (defined $remove_x_chr && length $remove_x_chr) {
            if ($remove_x_chr =~ /^(?:0|false|no|n)$/i) {
                push @cmd, '--no-remove-X-chr';
            }
            else {
                push @cmd, '--remove-X-chr';
            }
        }
        push @cmd, ('--display-gwas', $display_gwas) if length $display_gwas;
        push @cmd, ('--target-snps', $target_snps) if length $target_snps;
        push @cmd, ('--target-snp-genes', $target_snp_genes) if length $target_snp_genes;
        if (length $get_common_associations) {
            if ($get_common_associations =~ /^(?:1|true|yes|y)$/i) {
                push @cmd, '--get-common-associations';
            }
            else {
                push @cmd, "--get-common-associations=$get_common_associations";
            }
        }
        push @cmd, ('--common-association-top-hit-threshold', $common_association_top_hit_threshold)
          if length $common_association_top_hit_threshold;
        push @cmd, ('--reference-build', $reference_build) if length $reference_build;
        push @cmd, ('--local-gtf-window-bp', $local_gtf_window_bp) if length $local_gtf_window_bp;
        push @cmd, ('--local-max-hits-per-fig', $local_max_hits_per_fig)
          if defined $local_max_hits_per_fig && $local_max_hits_per_fig =~ /^\d+$/;
        push @cmd, ('--local-manhattan-columns', $local_manhattan_columns)
          if defined $local_manhattan_columns && $local_manhattan_columns =~ /^\d+$/;
        push @cmd, ('--local-manhattan-annotation', $local_manhattan_annotation)
          if length $local_manhattan_annotation;

        my $child_pid = fork();
        return { content => [{ type => 'text', text => "ERROR: fork failed" }] } unless defined $child_pid;

        if ($child_pid == 0) {
            open(my $child_out, '>', $out_file) or die "Failed to open $out_file: $!";
            open(STDIN, '<', File::Spec->devnull()) or die "Failed to redirect STDIN: $!";
            open(STDOUT, '>&', $child_out) or die "Failed to redirect STDOUT: $!";
            open(STDERR, '>&', \*STDOUT) or die "Failed to redirect STDERR: $!";
            exec @cmd;
            exit(1);
        }

        open(my $pfh, '>', $pid_file) or die "Failed to write pid file: $!";
        print $pfh "$child_pid\t$out_file\n";
        close $pfh;

        return {
            content => [{
                type => 'text',
                text => "QUERYING: gunplot wrapper started\nPID: $child_pid\nOutput file: $out_file\nAsk the AI agent to check status with: {\"spec_file\": \"$spec\", \"pid\": $child_pid}\nHints: use display_gwas=\"ALL_FEMALE\" for a single-GWAS plot set, display_gwas=\"EUR,EUR_FEMALE,EUR_MALE\" for mixed track display, and target_snps=\"rs123\" for inquiry-SNP local plots. In Ubuntu Docker, prefer local_manhattan/local_gtf first for a fast validation and add manhattan only when you want the slower genomewide render."
            }]
        };
    }
);

# # -------------------------
# $server->tool(
#     name        => 'run_cmd',
#     description => 'Run any user command in cmd',
#     input_schema => {
#         type       => 'object',
#         properties => {
#             user_cmd => { type => 'string' }
#         },
#         required => ['user_cmd']
#     },
#     code => sub ($tool, $args) {
#         my $cmd = $args->{user_cmd};
#         open(my $fh, '-|', $cmd) or return "Cannot execute command: $!";
#         local $/;
#         my $content = <$fh>;
#         close $fh;
#         return $content;
#     }
# );


# -------------------------
# get GTEx eQTL information for a given genesymbol, not working for gemini due to timeout but working for mcphost;
# -------------------------
# $server->tool(
#     name        => 'gtex_query_eQTLs4genesymbol',
#     description => 'Query GTEx eQTLs for a given gene symbol using the GTEx API, and return the results in a structured format. 
#     Note: this tool requires a local Perl script "gtex_query_eQTLs4genesymbol.pl" that handles the actual querying logic. 
#     The input is a gene symbol, and the output is the eQTL information related to that gene symbol, which can be used for downstream analysis
#     by other tools, such as the SAS tool defined above to run SAS codes via ODA on the queried eQTL information.',
#     input_schema => {
#         type       => 'object',
#         properties => {
#             genesymbol => { type => 'string', description => 'The gene symbol to query eQTLs for (e.g., "CD55")' },
#             output_file => { type => 'string', description => 'path to save the eQTL results in TSV format' }
#         },
#         required => ['genesymbol']
#     },
#     code => sub ($tool, $args) {
#         my $cmd = "perl -S gtex_query_eQTLs4genesymbol.pl --gene ".$args->{genesymbol};
#         if ($args->{output_file}) {
#             $cmd .= ">$args->{output_file}";

#         }
#         #print STDERR "Running GTEx eQTL query for gene symbol: $cmd\n";
#         ##Do not combine "-|" with the input $cmd directly, as it may cause issues with complex commands or arguments. Instead, use a safer approach to execute the command and capture its output.;
#         #open(my $fh, "-|", $cmd) or return "Cannot execute command: $!";
#         #local $/;
#         #my $content = <$fh>;
#         #close $fh;
#         #return $content;
#         #Avoiding outputting the content directly, as it may cause issues with MCP response parsing. Instead, we can return a success message and the path to the output file if provided.;
#         eval {
#             system($cmd) == 0 or die "Command failed with exit code: " . ($? >> 8);
#         };
#         return $@ ? "Error executing GTEx query: $@" : "GTEx eQTL query executed successfully for gene symbol: $args->{genesymbol}. Output saved to: " . ($args->{output_file} // "STDOUT");
#     }
# );

$server->tool(
    name        => 'gtex_async_query',
    description => 'Query GTEx for all eQTLs of a query gene by submitting the query as a background job, and return the results in a TSV file.'.
    'This tool is designed to avoid timeout issues when running long queries on platforms like Gemini. The tool checks the status of the '.
    'background job and returns partial output if still running, or the final results when complete.',
    input_schema => {
        type => 'object',
        properties => {
            query_genesymbol     => { type => 'string', 
                               description => 'Gene symbol (e.g., "CD55")' },
            output_file   => { type => 'string',
                               description => 'Optional output TSV file path (optional)' },
            pid           => { type => 'integer',
                               description => 'Optional PID to check status of a previous query (when run the tool again with the same query_snp)' }
        },
        required => ['query_genesymbol']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $query_genesymbol = $args->{query_genesymbol};
        #This output filename if provided with the previous pid, the tool will check the status of the previous query
        #and return the results in this file when the job is finished. If not provided, the tool will save the results 
        #in a default file named "tmp{timestamp}/haploreg_{query_snp}.tsv" as demonstrated in the later part of the code;
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        my $out_file; #Note: it is necessary to define the $out_file variable here to avoid of conflicts between global var and local var;
        my $pid_file; #Note: it is necessary to define the $pid_file variable here to avoid of conflicts between global var and local var;
        # If checking a previous PID
        if (defined $pid_arg) {
            $pid_file = "tmp*/GTEx_$query_genesymbol.pid";
            my @pid_files = glob($pid_file);
            
            foreach my $pf (@pid_files) {
                if (-f $pf) {
                    open(my $pfh, '<', $pf);
                    my $stored_pid = <$pfh>;
                    close $pfh;
                    chomp $stored_pid;
                    
                    if ($stored_pid == $pid_arg) {
                        #Note: the output file without "my" here is named "tmp{timestamp}/GTEx_{query_genesymbol}.tsv" by default 
                        #when the query is first run without the output_file argument, and the output file path 
                        #is saved in the same directory as the pid file. So when checking the status with the pid, 
                        #we can get the output file path by replacing the ".pid" extension
                        $out_file = $pf;
                        $out_file =~ s/\.pid$/.tsv/;
                        
                        my @pids = split(" ", `perl -S GetPIDs.pl`);
                        chomp @pids;
                        
                        if (grep { $_ == $pid_arg } @pids) {
                            return {
                                content => [{
                                    type => "text",
                                    text => "STATUS: RUNNING (PID $pid_arg)\n" .
                                           "Output file: $out_file\nAsk the AI agent to check status again in a moment."
                                }]
                            };
                        }
                        
                        # Job finished
                        my $content = "";
                        if (-f $out_file && open(my $fh, '<', $out_file)) {
                            local $/;
                            $content = <$fh>;
                            close $fh;
                        }
                        
                        unlink $pf;
                        cleanup_generated_tmpdir_if_empty($pf);
                        
                        return {
                            content => [{
                                type => "text",
                                text => "STATUS: COMPLETE (PID $pid_arg)\n\nResults saved to: $out_file\n\n" . $content
                            }]
                        };
                    }
                }
            }
            
            return {
                content => [{
                    type => "text",
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        # Start new query by creating a tmp directory and running the query in the background
        my $tmpdir = "tmp" . time();
        unless (-d $tmpdir) {
            mkdir $tmpdir or die "Failed to create tmp directory: $!";      
        }
        # Create a file named "tmp{timestamp}/GTEx_{query_genesymbol}.tsv" for first time query, 
        #and save the results in this file when the job is finished. If the tool is called again 
        #with the same query_genesymbol and the pid of the previous job, it will check the status of the 
        #previous job and return the results in this file when the job is finished;
        $pid_file = "$tmpdir/GTEx_$query_genesymbol.pid";

        #Note: the output file without "my" is named "tmp{timestamp}/GTEx_{query_genesymbol}.tsv" by default 
        #when the query is first run without the output_file argument,
        $out_file = $output_file // "$tmpdir/GTEx_$query_genesymbol.tsv";

        unless (-f $pid_file) {
            my $pid = fork();
            return { content => [{ type => "text", text => "ERROR: Could not fork." }] } unless defined $pid;

            if ($pid == 0) {
                my $cmd = "perl -S gtex_query_eQTLs4genesymbol.pl --gene '$query_genesymbol' "
                        . "> '$out_file'";
                print STDERR "Executing command in background: $cmd\n";
                exec "$cmd 2>&1";
                exit(1);
            }

            open(my $pfh, '>', $pid_file);
            print $pfh $pid;
            close $pfh;
            # Return initial status with PID and instructions for checking later, avoiding the endless loop of checking by the AI agent, 
            # which may cause excessive checking and potential issues with the server;
            return {
                content => [{
                    type => "text",
                    text => "QUERYING: GTEx for $query_genesymbol\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"query_genesymbol\": \"$query_genesymbol\", \"pid\": $pid}\n" . 
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking and only print the first 10 lines of output.\n"
                }]
            };
        }
    }
);

# $server->tool(
#     name        => 'gtex_async_query',
#     description => 'Query all eQTLs for a gene in GTEx asynchronously.',
#     input_schema => {
#         type => 'object',
#         properties => {
#             genesymbol  => { type => 'string' },
#             output_file => { type => 'string', 
#                              description => 'Output text file path for GTEx query results' }
#         },
#         required => ['genesymbol', 'output_file']
#     },
#     code => sub ($tool, $args) {

#         my $gene = $args->{genesymbol};
#         my $out  = $args->{output_file};
#         my $preview_lines = 100;
#         # Sanitize
#         $gene =~ s/[^a-zA-Z0-9_-]//g;
#         $out  =~ s/[^a-zA-Z0-9_\.-]//g;

#         my $pid_file = "$out.pid";

#         # --------------------------------------------------
#         # CASE 1: Job never started → START IT
#         # --------------------------------------------------
#         unless (-f $pid_file) {

#             my $pid = fork();
#             unless (defined $pid) {
#                 return {
#                     content => [{
#                         type => "text",
#                         text => "ERROR: Could not fork process."
#                     }]
#                 };
#             }

# if ($pid == 0) {

#     exec "perl -S gtex_query_eQTLs4genesymbol.pl --gene $gene > $out 2>&1";

#     exit(1);
# }

#             # --- PARENT ---
#             open(my $pfh, '>', $pid_file);
#             print $pfh $pid;
#             close $pfh;

#             return {
#                 content => [{
#                     type => "text",
#                     text => "STATUS: STARTED\nPID: $pid\nAsk AI agent to check the PID $pid later to get results."
#                 }]
#             };
#         }

#         # --------------------------------------------------
#         # CASE 2: Job exists → Check status
#         # --------------------------------------------------

#         open(my $pfh, '<', $pid_file);
#         my $pid = <$pfh>;
#         close $pfh;
#         chomp $pid;

#         # If still running
#         #my @pids=`ps -e -o pid=`;
#         my @pids=split(" ",`perl -S GetPIDs.pl`);
#         chomp @pids;
#         if (grep { $_ == $pid } @pids) {
#             my $preview = "";
#             if (-f $out && open(my $fh, '<', $out)) {
#                 local $/;
#                 $preview = substr(<$fh>, 0, $preview_lines);
#                 close $fh;
#             }

#             return {
#                 content => [{
#                     type => "text",
#                     text => "STATUS: RUNNING (PID $pid)\n\nPartial output:\n$preview\n\nCall again to check."
#                 }]
#             };
#         }

#         # --------------------------------------------------
#         # CASE 3: Job finished → Return result
#         # --------------------------------------------------

#         my $content = "";
#         if (-f $out && open(my $fh, '<', $out)) {
#             local $/;
#             $content = <$fh>;
#             close $fh;
#         }

#         unlink $pid_file;

#         return {
#             content => [{
#                 type => "text",
#                 text => "STATUS: COMPLETE\n\nResults are saved in $out\n\n"."Results (first $preview_lines lines):\n"
#                         . substr($content, 0, $preview_lines) 
#             }]
#         };
#     }
# );

#Query high LD SNPs in Haploreg4 and save the html file, which would fail when used by gemini due to timeout;
#But mcphost will not have the timeout issue when running this tool, which is kept here as a referene;
# $server->tool(
#     name        => 'haploreg4_fetch_and_parse_high_LD_snps',
#     description => 'Get high LD SNPs by quering HaploReg v4.1 and parses the stream directly into a TSV file.',
#     input_schema => {
#         type       => 'object',
#         properties => {
#             query_snp     => { type => 'string' },
#             ld_threshold  => { type => 'number'},
#             ld_population => { type => 'string', enum => ['AFR', 'AMR', 'ASN', 'EUR'] },
#             output_tsv    => { type => 'string' }
#         },
#         required => ['output_tsv']
#     },
#     code => sub ($tool, $args) {
#         my $query = $args->{query_snp} or die "No query_snp provided.";
#         my $output_tsv = $args->{output_tsv};
#         my $threshold = $args->{ld_threshold}  // 0.8;
#         my $pop       = $args->{ld_population} // 'EUR';

#         # Construct URL
#         my $url = sprintf(
#             "https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php?query=%s&ldThresh=%s&ldPop=%s",
#             $query, $threshold, $pop
#         );

#         # Pipe curl output directly into a filehandle
#         # --globoff handles special chars, -s is silent
#         open(my $pipe, "-|", "curl -s --globoff \"$url\"") or die "Could not open pipe to curl: $!";
        
#         # Slurp the streamed HTML
#         my $html = do { local $/; <$pipe> };
#         close($pipe);
#         return "Failed to fetch data from HaploReg4." unless $html;

#         # 1. Extract the Results Table
#         my ($table) = $html =~ /(<table[^>]*?resulttable.*?>.*?<\/table>)/is;
        
#         if (!$table) {
#             # Check for error messages in the stream
#             my ($error_msg) = $html =~ /<font color=red>(.*?)<\/font>/i;
#             die "HaploReg Error: " . ($error_msg // "Table not found. Check if the SNP exists.");
#         }

#         # 2. Extract Rows
#         my @rows = $table =~ /<tr.*?>(.*?)<\/tr>/gis;
#         my @header = ("chr", "pos_hg38", "r2", "D_prime", "variant", "Ref", "Alt", "AFR_freq", "AMR_freq", "ASN_freq", "EUR_freq", "Promoter", "Enhancer", "DNAse", "Proteins", "Motifs", "GWAS", "GRASP", "eQTL", "dbSNP_func");

#         # 3. Write to Output File
#         open(my $fh_out, '>', $output_tsv) or die "Cannot write to $output_tsv: $!";
#         print $fh_out join("\t", @header) . "\n";

#         my $count = 0;
#         foreach my $row (@rows) {
#             next if $row =~ /<b>chr<\/b>/i; # Skip the HTML header
#             my @cells = $row =~ /<td.*?>(.*?)<\/td>/gis;
#             next unless @cells;

#             my @clean_row = map {
#                 my $val = $_ // "";
#                 # rsID cleanup
#                 if ($val =~ />(rs\d+)</) {
#                     $val = $1;
#                 } else {
#                     $val = $1 if $val =~ /title="(.*?)"/i;
#                     $val =~ s/<.*?>//g;      # Strip tags
#                     $val =~ s/&nbsp;/ /g;    # Clean spaces
#                     $val =~ s/\s+/ /g;       # Collapse whitespace
#                     $val =~ s/^\s+|\s+$//g;  # Trim
#                 }
#                 $val;
#             } @cells;
            
#             # Ensure column count consistency
#             push @clean_row, ("") x (scalar(@header) - scalar(@clean_row)) if @clean_row < @header;
#             print $fh_out join("\t", @clean_row) . "\n";
#             $count++;
#         }
#         close($fh_out);

#         return "Successfully parsed $count rows into $output_tsv.";
#     }
# );

#Implementation of the above tool with asynchronous execution to avoid timeout issue, which can be used by gemini and mcphost;
#Use this as a template for other simplar tools that run system commands in background to avoid timeout issue when running on platforms like Gemini;
#promp for the AI agent to change the mcp tool for other similar tools that run system commands in background to avoid timeout issue 
#when running on platforms like Gemini is as follows:
#Based on the mcp tool provided at the beginning of these codes, please let it specifically as a tool for the perl script, the contents for query_high_LD_SNPs_at_Haploreg4.pl is provided at the bottom.

$server->tool(
    name        => 'query_haploreg',
    description => 'Query HaploReg4 for high LD SNPs by submitting the query as a background job, and return the results in a TSV file.'.
    'This tool is designed to avoid timeout issues when running long queries on platforms like Gemini. The tool checks the status of ' .
    'the background job and returns partial output if still running, or the final results when complete.',
    input_schema => {
        type => 'object',
        properties => {
            query_snp     => { type => 'string', 
                               description => 'SNP identifier (e.g., rs17425819)' },
            ld_threshold  => { type => 'number',
                               description => 'LD threshold (default: 0.8)',
                               default => 0.8 },
            ld_population => { type => 'string',
                               description => 'Population: AFR, AMR, ASN, EUR (default: EUR)',
                               default => 'EUR',
                               enum => ['AFR', 'AMR', 'ASN', 'EUR'] },
            output_file   => { type => 'string',
                               description => 'Optional output TSV file path (optional)' },
            pid           => { type => 'integer',
                               description => 'Optional PID to check status of a previous query (when run the tool again with the same query_snp)' }
        },
        required => ['query_snp']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $query_snp = $args->{query_snp};
        my $ld_threshold = $args->{ld_threshold} // 0.8;
        my $ld_population = $args->{ld_population} // 'EUR';
        #This output filename if provided with the previous pid, the tool will check the status of the previous query
        #and return the results in this file when the job is finished. If not provided, the tool will save the results 
        #in a default file named "tmp{timestamp}/haploreg_{query_snp}.tsv" as demonstrated in the later part of the code;
        my $output_file = $args->{output_file};
        my $pid_arg = $args->{pid};
        my $out_file; #Note: it is necessary to define the $out_file variable here to avoid 
        #the "Use of uninitialized value $out_file in concatenation (.) or string at line codes" 
        #error when the tool is called for the first time without the pid argument,
        # If checking a previous PID
        my $pid_file; #Note: it is necessary to define the $pid_file variable here to avoid of conflicts between global var and local var;
        if (defined $pid_arg) {
            $pid_file = "tmp*/haploreg_$query_snp.pid";
            my @pid_files = glob($pid_file);
            
            foreach my $pf (@pid_files) {
                if (-f $pf) {
                    open(my $pfh, '<', $pf);
                    my $stored_pid = <$pfh>;
                    close $pfh;
                    chomp $stored_pid;
                    
                    if ($stored_pid == $pid_arg) {
                        #Note: when checking the status with the pid, we can get the output file path
                        #by replacing the ".pid" extension in the pid file name with ".tsv".
                        #This way, we can ensure that we are checking the correct output file associated with the specific PID.
                        $out_file = $pf;
                        $out_file =~ s/\.pid$/.tsv/;
                        
                        my @pids = split(" ", `perl -S GetPIDs.pl`);
                        chomp @pids;
                        
                        if (grep { $_ == $pid_arg } @pids) {
                            return {
                                content => [{
                                    type => "text",
                                    text => "STATUS: RUNNING (PID $pid_arg)\n" .
                                           "Output file: $out_file\nAsk the AI agent to check status again in a moment."
                                }]
                            };
                        }
                        
                        # Job finished
                        my $content = "";
                        if (-f $out_file && open(my $fh, '<', $out_file)) {
                            local $/;
                            $content = <$fh>;
                            close $fh;
                        }
                        
                        unlink $pf;
                        cleanup_generated_tmpdir_if_empty($pf);
                        
                        return {
                            content => [{
                                type => "text",
                                text => "STATUS: COMPLETE (PID $pid_arg)\n\nResults saved to: $out_file\n\n" . $content
                            }]
                        };
                    }
                }
            }
            
            return {
                content => [{
                    type => "text",
                    text => "ERROR: PID $pid_arg not found or already completed."
                }]
            };
        }

        # Start new query by creating a tmp directory and running the query in the background
        my $tmpdir = "tmp" . time();
        unless (-d $tmpdir) {
            mkdir $tmpdir or die "Failed to create tmp directory: $!";      
        }
        # Create a file named "tmp{timestamp}/haploreg_{query_snp}.tsv" for first time query, 
        #and save the results in this file when the job is finished. If the tool is called again 
        #with the same query_snp and the pid of the previous job, it will check the status of the 
        #previous job and return the results in this file when the job is finished;
        $pid_file = "$tmpdir/haploreg_$query_snp.pid";
        $out_file = $output_file // "$tmpdir/haploreg_$query_snp.tsv";

        unless (-f $pid_file) {
            my $pid = fork();
            return { content => [{ type => "text", text => "ERROR: Could not fork." }] } unless defined $pid;

            if ($pid == 0) {
                my $cmd = "perl -S query_high_LD_SNPs_at_Haploreg4.pl "
                        . "--query-snp '$query_snp' "
                        . "--ld-threshold $ld_threshold "
                        . "--ld-population $ld_population "
                        . "--output-tsv '$out_file'";
                
                exec "$cmd 2>&1";
                exit(1);
            }

            open(my $pfh, '>', $pid_file);
            print $pfh $pid;
            close $pfh;
            # Return initial status with PID and instructions for checking later, avoiding the endless loop of checking by the AI agent, 
            # which may cause excessive checking and potential issues with the server;
            return {
                content => [{
                    type => "text",
                    text => "QUERYING: HaploReg4 for $query_snp\nPID: $pid\nOutput file: $out_file\n" .
                           "Ask the AI agent to check status with: {\"query_snp\": \"$query_snp\", \"pid\": $pid}\n" . 
"The loop should be no more than 3 times with an interval of 30 seconds to avoid excessive checking and only print the first 10 lines of output.\n"
                }]
            };
        }
    }
);


##########################################################System tools for debugging and environment checks################################################

# -------------------------
# OS information
# -------------------------
$server->tool(
    name        => 'uname',
    description => 'Get OS and kernel information',
    input_schema => {
        type       => 'object',
        properties => {}
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        return `uname -a`;
    }
);

# -------------------------
# Environment variables
# -------------------------
$server->tool(
    name        => 'env_vars',
    description => 'List environment variables',
    input_schema => {
        type       => 'object',
        properties => {}
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        return join("\n", map { "$_=$ENV{$_}" } sort keys %ENV);
    }
);

# -------------------------
# Disk usage
# -------------------------
$server->tool(
    name        => 'disk_usage',
    description => 'Show disk usage (df -h)',
    input_schema => {
        type       => 'object',
        properties => {}
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        return `df -h`;
    }
);

# -------------------------
# Uptime
# -------------------------
$server->tool(
    name        => 'uptime',
    description => 'Show system uptime',
    input_schema => {
        type       => 'object',
        properties => {}
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        return `uptime`;
    }
);

#An internal helper function to write content into a file, which can be used by other tools defined in this server;
sub internal_write_file {
    my ($path, $content) = @_;
    open(my $fh, '>', $path) or return "Error: $!";
    print $fh $content;
    close $fh;
    return "Successfully wrote to $path";
}


# Example usage
# Note: you need to have the 'pull_column' script in your PATH or current directory for this to work.
# my $is_perl_script_runnable = is_perl_script_runnable("pull_column");
# print "Perl script is runnable: " . ($is_perl_script_runnable ? "Yes" : "No") . "\n";
sub is_perl_script_runnable {
    my ($script_name) = @_;
    
    # Use perl -S to check if the script is runnable
    my $command = "perl -S $script_name";
    my $output = `$command 2>&1`;
    my $exit_code = $? >> 8;
    
    # Return true if the script is runnable (exit code 0)
    return $exit_code == 0;
}

#Check whether the specific perl script is runnable;
$server->tool(
    name        => 'is_perl_script_runnable',
    description => 'Check if a Perl script is runnable using perl -S',
    input_schema => {
        type       => 'object',
        properties => {
            script_name => { type => 'string', description => 'Name of the Perl script to check' }
        },
        required => ['script_name']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $script_name = $args->{script_name};
        return is_perl_script_runnable($script_name) ? "Yes" : "No";
    }
);

$server->tool(
    name        => 'check_perl_script_functionality',
    description => 'Check the functionality of a Perl script by attempting to run it with perl -S',
    input_schema => {
        type       => 'object',
        properties => {
            script_name => { type => 'string', description => 'Name of the Perl script to check' }
        },
        required => ['script_name']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $script_name = $args->{script_name};
        my $command = "perl -S $script_name"; # Attempt to run the perl script
        my $output = `$command 2>&1`;
        return $output ? $output : "Script ran successfully with no output.";   
    }
);

#This is a sample tool definition for an echo tool.
# It is commented out by default, because it will interupt other tools' functionality if left enabled.
# for example, it will affect the 'write_file' tool's ability to return structured MCP responses.
#as the mcphost would run this echo tool for every input, causing unexpected behavior.
# ------------------------- # Echo # ------------------------- 
#$server->tool( 
#name => 'echo', 
#description => 'Echo the input text', 
#input_schema => { 
#type => 'object', 
#properties => { msg => { type => 'string' } }, 
#required => ['msg'] }, 
#code => sub ($tool, $args) {
# return "Echo: $args->{msg}"; 
#} 
#);

# -------------------------
# System time
# -------------------------
$server->tool(
    name        => 'system_time',
    description => 'Get current system time',
    input_schema => {
        type       => 'object',
        properties => {}
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        return strftime('%Y-%m-%d %H:%M:%S %Z', localtime);
    }
);

# -------------------------
# Current working directory
# -------------------------
$server->tool(
    name        => 'pwd',
    description => 'Get current working directory',
    input_schema => {
        type       => 'object',
        properties => {}
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        return getcwd();
    }
);

# -------------------------
# List directory contents
# -------------------------
$server->tool(
    name        => 'list_dir',
    description => 'List contents of a directory',
    input_schema => {
        type       => 'object',
        properties => {
            path => { type => 'string', description => 'The directory path to list' }
        }
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $path = $args->{path} // '.';
        opendir(my $dh, $path) or return "Cannot open directory: $!";
        my @files = readdir($dh);
        closedir $dh;
        return join("\n", sort @files);
    }
);

# -------------------------
# Read text file
# -------------------------
$server->tool(
    name        => 'read_file',
    description => 'Read a text or tsv file',
    input_schema => {
        type       => 'object',
        properties => {
            path => { type => 'string' }
        },
        required => ['path']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $path = $args->{path};
        open(my $fh, '<', $path) or return "Cannot open file: $!";
        local $/;
        my $content = <$fh>;
        close $fh;
        return $content;
    }
);

# -------------------------
# write text file
# -------------------------
$server->tool(
    name        => 'write_file',
    description => 'Write text content to a file',
    input_schema => {
        type       => 'object',
        properties => {
            outpath => {
                type        => 'string',
                description => 'Output file path (relative or absolute)'
            },
            content => {
                type        => 'string',
                description => 'Text content to write into the file'
            },
            overwrite => {
                type        => 'boolean',
                description => 'Whether to overwrite the file if it exists'
            }
        },
        required => ['outpath', 'content']
    },
    code => sub ($tool, $args) {
        cleanup_all_generated_empty_tmpdirs();
        my $path      = $args->{outpath};
        my $content   = $args->{content};
        my $overwrite = $args->{overwrite} // 1;

        # 1. Check if file exists and overwrite is forbidden
        if (-e $path && !$overwrite) {
            return {
                isError => \1, # JSON boolean true
                content => [{ type => 'text', text => "Error: File already exists at $path" }]
            };
        }

        # 2. Attempt to write the file
        # Added ':utf8' to ensure correct encoding
        #if replacing >:utf8 with >, it may cause issues with non-ASCII characters;
        if (open(my $out_fh, '>', $path)) {
            print {$out_fh} $content;
            close $out_fh;

            # 3. Return the standard MCP Result structure
            return {
                content => [{
                    type => 'text',
                    text => "Successfully wrote " . length($content) . " bytes to $path"
                }]
            };
        } else {
            # 4. Handle OS-level write errors (permissions, path not found, etc.)
            return {
                isError => \1,
                content => [{ type => 'text', text => "System Error: $!" }]
            };
        }
    }
);

# # -------------------------
# # fetch, clean, and save web content, which has been implemented in the perl script: fetch_a_website;
# # which can be run it directly in the command line;
# # -------------------------
# $server->tool(
#     name        => 'fetch_clean_save_web',
#     description => 'Fetches a URL, cleans tracking parameters, strips heavy tags, and saves to file.',
#     input_schema => {
#         type       => 'object',
#         properties => {
#             url => { 
#                 type        => 'string', 
#                 description => 'The URL to fetch' 
#             },
#             output_file => { 
#                 type        => 'string', 
#                 description => 'Path to save the cleaned HTML' 
#             },
#             strip_scripts => { 
#                 type    => 'boolean', 
#                 description => 'If true, removes <script> and <style> tags' 
#             }
#         },
#         required => ['url', 'output_file'],
#     },
#     code => sub ($tool, $args) {
#         my $url           = $args->{url};
#         my $output_file   = $args->{output_file};
#         my $strip_scripts = $args->{strip_scripts} // 1;

#         # --- 1. Clean the URL (Removing tracking junk) ---
#         # Removes utm_*, ref, source, etc. from the query string
#         $url =~ s/([?&])(utm_[^&]+|ref=[^&]+|source=[^&]+|fbclid=[^&]+)(&?)/$3 ? $1 : ""/ge;
#         $url =~ s/\?$//; # Remove trailing question mark if query is now empty

#         # --- 2. Fetch the Content ---
#         my $http = HTTP::Tiny->new(
#             agent      => 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)',
#             verify_SSL => 0 # Recommended for Windows to avoid cert store issues
#         );
#         my $response = $http->get($url);

#         if ($response->{success}) {
#             my $content = $response->{content};

#             # --- 3. Clean the HTML ---
#             if ($strip_scripts) {
#                 # Remove <script>...</script> blocks
#                 $content =~ s/<script\b[^>]*>.*?<\/script>//gis;
#                 # Remove <style>...</style> blocks
#                 $content =~ s/<style\b[^>]*>.*?<\/style>//gis;
#                 # Remove HTML comments
#                 $content =~ s///gs;
#             }

#             # --- 4. Save to File ---
#             if (open(my $fh, '>', $output_file)) {
#                 binmode($fh, ":utf8"); # Crucial for Windows encoding
#                 print $fh $content;
#                 close($fh);
                
#                 return "Successfully Saved cleaned HTML to: $output_file";
#             } else {
#                 die "Could not open '$output_file': $!";
#             }
#         } else {
#             die "Fetch failed: " . ($response->{reason} // "Unknown error");
#         }
#     }
# );


# -------------------------
# MCP endpoint
# -------------------------
#any '/mcp' => $server->to_action;
# ---------------------------------------------------------
# ---------------------------------------------------------
# 3. GLOBAL SESSION INJECTOR HOOK
# ---------------------------------------------------------
app->hook(before_dispatch => sub ($c) {
    my $json = $c->req->json;
    return unless $json;
    if (my $sid = $c->param('session_id')) {
        $json->{params} //= {};
        $json->{params}{sessionId} = $sid;
        $c->req->json($json);
    }
});

# ---------------------------------------------------------
# 4. MAIN ROUTE
# ---------------------------------------------------------
any '/mcp' => sub ($c) {
    my $json = $c->req->json;
    return $c->render(json => { error => "No JSON" }, status => 400) unless $json;

    my $method = $json->{method} // '';

    # Log every incoming request for debugging    
    #warn "[MCP] Method: $method | Raw: " . (eval { encode_json($json) } // 'undef') . "\n";
    # --------------------------------------------------
    # 4a. INITIALIZE — generate & track a real session
    # --------------------------------------------------
    if ($method eq 'initialize') {
        my $sid = md5_sum(time() . rand() . $$);
        $sessions{$sid} = { created => time() };
        warn "[MCP] New session: $sid\n";

        return $c->render(json => {
            jsonrpc => "2.0",
            id      => $json->{id},
            result  => {
                protocolVersion => "2024-11-05",
                capabilities    => { tools => { listChanged => \1 } },
                serverInfo      => { name => "PerlBioServer", version => "1.0.0" },
                sessionId       => $sid
            }
        });
    }

    # --------------------------------------------------
    # 4b. notifications/initialized — handshake ACK
    # --------------------------------------------------
    if ($method eq 'notifications/initialized') {
        #warn "[MCP] Received initialized notification\n";
        return $c->render(json => { jsonrpc => "2.0", result => {} });
    }

    # --------------------------------------------------
    # 4c. TOOLS/LIST — filter out null/invalid entries
    #     MCP::Server stores internal metadata in the
    #     tools array; skip anything missing required fields.
    # --------------------------------------------------
    if ($method eq 'tools/list') {
        my @tool_list;
        if (my $tools_array = $server->{tools}) {
            foreach my $tool (@$tools_array) {
                #print "\n", Dumper($tool),"\n" unless (defined $tool->{name} || $tool->{description});
                next unless (defined $tool->{name} || $tool->{description});

                my $schema = $tool->{input_schema} || $tool->{inputSchema};
                # next unless ref($schema) eq 'HASH';

                push @tool_list, {
                    name        => $tool->{name},
                    description => $tool->{description},
                    inputSchema => $schema
                };
            }
        }

        warn "[MCP] Returning " . scalar(@tool_list) . " tools\n";

        return $c->render(json => {
            jsonrpc => "2.0",
            id      => $json->{id},
            result  => { tools => \@tool_list }
        });
    }

    # --------------------------------------------------
    # 4d. TOOLS/CALL — execute directly, no session check
    #     mcphost/gemini do NOT resend session ID on calls.
    # --------------------------------------------------
    if ($method eq 'tools/call') {

        my $tool_name = $json->{params}{name}      // '';
        my $tool_args = $json->{params}{arguments} // {};

        #warn "[MCP] Tool call: $tool_name\n";

        # Find the tool — use same filtering as tools/list
        my ($tool) = grep {
            defined $_->{name}
            && $_->{name} eq $tool_name
        } @{ $server->{tools} };

        unless ($tool) {
            return $c->render(json => {
                jsonrpc => "2.0",
                id      => $json->{id},
                error   => { code => -32601, message => "Unknown tool: $tool_name" }
            });
        }

        my $result = eval { $tool->{code}->($tool, $tool_args) };
        if ($@) {
            return $c->render(json => {
                jsonrpc => "2.0",
                id      => $json->{id},
                error   => { code => -32603, message => "Execution error: $@" }
            });
        }

        my $text = ref($result) ? encode_json($result) : ($result // '');

        return $c->render(json => {
            jsonrpc => "2.0",
            id      => $json->{id},
            result  => { content => [{ type => "text", text => $text }] }
        });
    }

    # --------------------------------------------------
    # 4e. CATCH-ALL — ACK notifications, error on unknown
    #     requests
    # --------------------------------------------------
    if (!defined $json->{id}) {
        warn "[MCP] Unhandled notification: $method\n";
        return $c->render(json => { jsonrpc => "2.0", result => {} });
    }

    warn "[MCP] Unhandled method: $method\n";
    return $c->render(json => {
        jsonrpc => "2.0",
        id      => $json->{id},
        error   => { code => -32601, message => "Method not found: $method" }
    });
    
    #No need to run this!
    # Fallback for anything else
    #return $server->to_action->($c);
};

#any '/mcp' => $server->to_action;

# ---------------------------------------------------------
# 5. START
# ---------------------------------------------------------
{
    local *STDOUT = $MCP_OUT;
    app->start;
}
