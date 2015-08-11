from ruffus import *
from seq_pipe import utils


@collate(input_files, formatter("([^/]+)_[12].fastq$"), ["{path[0]}/{1[0]}_1.fastq", "{path[0]}/{1[0]}_2.fastq"])
def collate_files(input_files, output_files):
    log.info("Collating paired fastq files: \n\t{} \n\t{}\n".format(input_files[0], input_files[1]))


@transform(collate_files, formatter("([^/]+)_[12].fastq$"), options.output+"{1[0]}.assembled.fastq", extras)
def pear_fastq_files(input_files, output_file, extras):

    log.info("Starting pear run on %s and %s", input_files[0], input_files[1])

    output_file = re.sub(r"\.assembled\.fastq", "", output_file)
    args = ["pear", "-f", input_files[0], "-r", input_files[1], "-o", output_file]

    utils.run_cmd("pear")

@transform(input_files, suffix(".fastq"), ".fastq", extras)
def upload_to_one_codex(input_file, output_file, extras):

    args = ["onecodex", "upload", input_file]
    log.info("uploading %s to One Codex", input_file)

    utils.run_cmd(args, "One Codex")

# only bowtie
# rsem does this for you
@transform(rename_accepted_hits, suffix(".bam"),".sorted.bam", extras)
def sort_bam(input_file, output_file, extras):
    log.info("Sorting %s ", input_file)

    # hacky
    output_file = re.sub(r"\.bam", "", output_file)
    args = ["samtools-rs", "rocksort", "-@", "8", "-m", "16G", input_file, output_file]):
    utils.run_cmd(args, "samtools rocksort")

    # careful
    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)


@transform(sort_bam, suffix(".sorted.bam"), ".bed", options.output, extras)
def bam_to_bed(input_file, output_file, output, extras):

    log.info("Converting %s to a bed file", input_file)
    args = ["bamToBed", "-i", input_file ">", output_file]
    utils.run_cmd(args, "bam_to_bed")

    # now we can move sorted bam to output
    file_name = os.path.basename(input_file)
    new_name = os.path.join(output, file_name)
    os.rename(input_file, new_name)


@transform(bam_to_bed, suffix(".bed"), ".bg", options.size, extras)
def bed_to_bg(input_file, output_file,  size_file, extras):
    log.info("Converting %s to a genome coverage file", input_file)
    
    args = ["genomeCoverageBed", "-bg", "-split", "-i", input_file, "-g", size_file, ">", output_file]
    utils.run_cmd(args, "bed_to_bg")

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)


@transform(bed_to_bg, suffix(".bg"), ".bw", genome, options.output)
def bg_to_bw(input_file, output_file, genome, output):
    log.info("Creating bigwig file from bg:  %s", input_file)
    command = "bedGraphToBigWig {} {} {}".format( input_file, size_file, output_filet)
    
    if subprocess.call(command, shell=True):
        log.warn("bg to bw conversion of %s failed, exiting", input_file)
        extras.report_error("bg_to_bw","bg to bw conversion of {} failed".format(input_file))
        raise SystemExit

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)


@transform(input_files, suffix(".fastq"), ".sam", options, stats_file, extras)
def align_with_bowtie_two(input_file, output_file, options, stats_file, extras):
    
    log.info("Running bowtie2 on %s", input_file)
    
    # use poen explicitly to cpature STDERR, check still
    args = ["bowtie2", "-t", "--no-unal", "-p", str(options.cores), "-x", options.index, input_file, "-S", output_file]
    output = utils.run_cmd(args, "bowtie2", True)

    # pass along to be saved
    utils.bowtie_record_output(output, input_file, stats_file)


# call out to external bwtools here
@merge(bg_to_bw, os.path.join(options.output,"bigWigStats-"+time_stamp+".out"))
def bw_stats(input_files, output_file):
        
    # we are going to call bwtool summary and bwtool distribution
    # have to explicitly send stdout stuff like that
    # what a program
    summary = "bwtool summary 10000 -header -with-sum {} /dev/stdout"
    dist    = "bwtool distribution {} /dev/stdout"
    

    for input_file in input_files:
        log.info("Running bigwig stats on {}".format(input_file))
        with open(output_file, "a+") as stats:
            for command in [summary, dist]:

                output = utils.run_cmd(command.format(os.path.abspath(input_file)).split(), command, True)

                if command.startswith("bwtool summary"):
                    stats.write("#### bwtool summary for {}\n".format(input_file))
                    stats.write(output)
                    stats.write("####\n")
          
                # filter zeros out
                else:
                    output = output.rstrip()
                    output_clean = [line for line in output.split("\n") if line.split('\t')[1] != '0']
                    stats.write("#### bwtool distribution for {}\n".format(input_file))
                    stats.write("depth\tcount\n")
                    stats.write("\n".join(output_clean))
                    stats.write("\n####\n")

            stats.write("\n\n")


@transform(input_files, formatter(), options.output+"{basename[0]}.genes.results", options, extras)
def rsem_align(input_file, output_file, options, extras):
    
    mean_len = extras.get_mean_length(input_file)
    output_file = output_file.replace(".genes.results", "")

    log.info("Running rsem calc exp on %s", output_file)
    command = ["rsem-calculate-expression", "-p", str(options.cores), "--calc-ci", input_file, options.index, output_file]

    run_cmd(command)


@transform(rsem_align, suffix(".genes.results"), ".pdf", options)
def rsem_plot_model(input_file, output_file, options):

    sample_name = input_file.replace("\.genes\.results", "")
    command = ["rsem-plot-model", sample_name, output_file]

    run_cmd(command)


@merge(rsem_align, "gene_exp.mtx", extras)
def rsem_generate_exp_matrix(input_files, matrix, extras):

    sample_list = extras.gen_sample_list()
    command = ["rsem-generate-data-matrix", sample_list, ">", matrix]

    run_cmd(command)


@transform(generate_exp_matrix, suffix(".mtx"), ".diff", extras)
def rsem_run_ebseq(input_file, output_file, extras):

    cond_str = extras.gen_cond_str()
    command = ["rsem-run-ebseq", input_file, cond_str, output_file]

    run_cmd(command)


@transform(run_ebseq, suffix(".diff"), ".sigdiff", options)
def rsem_fdr_correct(input_file, output_file, options):

    command = ["rsem-control-fdr", input_file, str(options.fdr), output_file]

    run_cmd(command)


@transform(write_excel_sheet, suffix(".xlsx"), ".email", options, input_files, extras)
def report_success(input_file, output_file, options, inputfiles, extras):
    
    log.info("Sending email report")
    
    # Create a text/plain message
    email_body = []
    email_body.append("Differential expression pipeline results:\n")
    email_body.append("The following fastq files were used:")
    
    for file in input_files:
        email_body.append("- {}".format(file))

    email_body.append("\nThe results (xlsx spreadsheet) and pipeline log are attatched")
    email_body.append("Please direct any questions to kgmcchesney@wisc.edu")

    # msg object
    msg = MIMEMultipart()

    # header stuff
    # no one else cares but me!
    root  = "root@alpha-helix.oncology.wisc.edu"
    subject = "Tophat DE pipeline Success report: {}".format(time.strftime("%d/%m/%Y"))

    msg['Subject'] = subject
    msg['From'] = root
    msg['To'] = COMMASPACE.join(options.emails)
    msg.attach( MIMEText("\n".join(email_body)) )

    # attatch the files
    for file in [input_file, log.handlers[0].baseFilename]:
        part = MIMEBase('application', "octet-stream")
        part.set_payload( open(file,"rb").read() )
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(file))
        msg.attach(part)
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(root, options.emails, msg.as_string())
    s.quit()
