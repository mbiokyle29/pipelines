from ruffus import *

@transform(input_files, suffix(".fastq"), ".fastq", extras)
def upload_to_one_codex(input_file, output_file, extras):

    args = ["onecodex", "upload", input_file]
    log.info("uploading %s to One Codex", input_file)

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("One Codex","uploading to One Codex failed \n{}".format(output))
        raise SystemExit


@collate(input_files, formatter("([^/]+)_[12].fastq$"), ["{path[0]}/{1[0]}_1.fastq", "{path[0]}/{1[0]}_2.fastq"])
def collate_files(input_files, output_files):
    log.info("Collating paired fastq files: \n\t{} \n\t{}\n".format(input_files[0], input_files[1]))

@transform(collate_files, formatter("([^/]+)_[12].fastq$"), options.output+"{1[0]}.assembled.fastq", extras)
def pear_fastq_files(input_files, output_file, extras):

    log.info("Starting pear run on %s and %s", input_files[0], input_files[1])

    output_file = re.sub(r"\.assembled\.fastq", "", output_file)
    args = ["pear", "-f", input_files[0], "-r", input_files[1], "-o", output_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("PEAR","Merging paried files with PEAR failed: \n{}".format(output))
        raise SystemExit

@transform(pear_fastq_files, suffix(".fastq"), ".fastq", extras)
def upload_paired_to_one_codex(input_file, output_file, extras):

    args = ["onecodex", "upload", input_file]
    log.info("uploading %s to One Codex", input_file)

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("One Codex","uploading peared files to One Codex failed \n{}".format(output))
        raise SystemExit


# paired alignment
@transform(collate_files, formatter("([^/]+)_[12].fastq$"), options.output+"{1[0]}-tophat-results/accepted_hits.bam", options, extras)
def tophat_align_paired(input_files, output_file, options, extras):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    args = ["tophat2", "-G", options.gtf,"-p", str(options.cores), "-o", output, options.index, input_files[0], input_files[1]]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("tophat_paired","tophat2 paired run failed:\n{}".format(output))
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

# unpaired alignment function
@transform(input_files, formatter(), options.output+"{basename[0]}-tophat-results/accepted_hits.bam", options, extras)
def tophat_align_unpaired(input_file, output_file, options, extras):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    log.info("Starting tophat2 run on %s", input_file)
    args = ["tophat2", "-G", options.gtf,"-p", options.cores, "-o", output, options.index, input_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn(output)
        extras.report_error("tophat_unpaird","tophat2 unpaired failed: \n{}".format(output))
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

# get rid of stupid accepted_hits.bam file
@transform([tophat_align_unpaired, tophat_align_paired],
            formatter(r".*/([^/]+)-tophat-results/accepted_hits.bam$"), 
            ''.join([options.output, "{1[0]}", ".bam"]), extras )
def rename_accepted_hits(input_file, output_file, extras):

    try:
        os.rename(input_file, output_file)
    except OSError:
        extras.report_error("rename_accepted_hits","Renaming {} to {} failed".format(input_file, output_file),log)

# both the tophat functions can feed into here
@transform(rename_accepted_hits, suffix(".bam"),".sorted.bam", extras)
def sort_bam(input_file, output_file, extras):
    log.info("Sorting %s ", input_file)

    # hacky
    output_file = re.sub(r"\.bam", "", output_file)

    if subprocess.call(["samtools-rs", "rocksort", "-@", "8", "-m", "16G", input_file, output_file]):
        log.warn("bam sorting %s failed, exiting", input_file)
        extras.report_error("sort_bam","bam sorting {} failed".format(input_file))
        raise SystemExit
    
    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

@transform(sort_bam, suffix(".sorted.bam"), ".bed", options.output, extras)
def bam_to_bed(input_file, output_file, output, extras):

    log.info("Converting %s to a bed file", input_file)
    if subprocess.call("bamToBed -i {} > {}".format(input_file, output_file), shell=True):
        log.warn("bam to bed conversion of %s failed, exiting", input_file)
        extras.report_error("bam_to_bed","bam to bed conversion of {} failed".format(input_file))
        raise SystemExit

    # now we can move sorted bam to output
    file_name = os.path.basename(input_file)
    new_name = os.path.join(output, file_name)
    os.rename(input_file, new_name)

@transform(bam_to_bed, suffix(".bed"), ".bg", options.size, extras)
def bed_to_bg(input_file, output_file,  size_file, extras):
    log.info("Converting %s to a genome coverage file", input_file)
    
    command = "genomeCoverageBed -bg -split -i {} -g {} > {}".format( input_file, size_file, output_file)
    if subprocess.call(command, shell=True):
        log.warn("bed to coverage conversion of %s failed, exiting", input_file)
        extras.report_error("bed_to_bg","bed to bg conversion of {} failed".format(input_file))
        raise SystemExit

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

@merge(sort_bam, os.path.join(options.output,"cuffdiff/gene_exp.diff",), options, extras)
def run_cuffdiff(input_files, output_file, options, extras):
    
    # grab the conf
    conf = extras.process_de_conf(options.de_conf)

    # make the output file
    output_dir = os.path.join(options.output,"cuffdiff/")

    # get the label
    neg_label = conf['negative-condition']['label']
    pos_label = conf['positive-condition']['label']
    label_string = "{},{}".format(neg_label, pos_label)

    neg_files = []
    pos_files = []

    # match the files in conf to the bams we have
    for conf_file in conf['negative-condition']['files']:

        conf_basename = os.path.splitext(conf_file)[0]
        # check the input files
        for bam_file in input_files:
            bam_basename = os.path.basename(bam_file).split(os.extsep)[0]
            
            # see if they match
            if conf_basename == bam_basename:
                log.info("%s matched %s in negative-condition", conf_file, bam_file)
                neg_files.append(bam_file)

                # remove it for optimization :) I think
                input_files.remove(bam_file)

     # match the files in conf to the bams we have
    for conf_file in conf['positive-condition']['files']:
        conf_basename = os.path.splitext(conf_file)[0]

        # check the input files
        for bam_file in input_files:
            bam_basename = os.path.basename(bam_file).split(os.extsep)[0]

            # see if they match
            if conf_basename == bam_basename:
                log.info("%s matched %s in positive-condition", conf_file, bam_file)
                pos_files.append(bam_file)

                # remove it for optimization :) I think
                input_files.remove(bam_file)


    # make sure both groups got something and that the inputfiles is emtpy
    if len(neg_files) == 0 or len(pos_files) == 0 or len(input_files) != 0:
        log.warn("Cuffdiff file specification error!")
        extras.report_error("cuffdiff", "Cuffdiff file specification error!")
        raise SystemExit

    # call it
    log.info("Starting cuffdiff run")
    args = ["cuffdiff", "-o", output_dir, "-L", label_string, "-p", options.cores, options.gtf, ','.join(neg_files), ','.join(pos_files)]
    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("Cuffdiff failed: %s", output)
        extras.report_error("cuffdiff","Cuffdiff failed: \n{}".format(output))
        raise SystemExit

    # print output
    log.info("cuffdiff output:")
    log.info(output)

@transform(run_cuffdiff, suffix(".diff"), ".xlsx", options, extras)
def write_excel_sheet(input_file, output_file, options, extras):
    
    # best message ever
    log.info("Writing gene expression results to spreadsheets")
    keep_cols = ["test_id", "locus", "sample_1", "sample_2", "status", "value_1",
                "value_2", "log2(fold_change)", "test_stat", "p_value", "q_value", "significant"]

    # read it in
    df = pd.read_table(input_file, sep="\t", usecols=keep_cols)

    # change the na,e of test_id
    keep_cols[0] = "gene_name"
    df.columns = keep_cols

    # annotate the genes
    if(options.annotation_db or options.annotation_file):
        df['gene_annotation'] = df['gene_name'].apply(extras.annotate_gene)

    # sort sig to the top
    df.sort(columns="significant", ascending=False, inplace=True, kind="heapsort")

    # write the sucker to a file
    writer = pd.ExcelWriter(output_file)
    df.to_excel(writer, sheet_name="Sheet1", index=False)
    writer.save()

    log.info("results written to: %s", output_file)

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

@transform(input_files, suffix(".fastq"), ".sam", options, stats_file, extras)
def align_with_bowtie(input_file, output_file, options, stats_file, extras):
    log.info("Running bowtie2 on %s", input_file)
    
    # use poen explicitly to cpature STDERR, check still
    args = ["bowtie2", "-t", "--no-unal", "-p", str(options.cores), "-x", options.index, input_file, "-S", output_file]
    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("Bowtie failed")
        raise SystemExit

    # print output
    log.info("Bowtie output:")
    log.info(output)

    # pass along to be saved
    extras.record_bowtie_output(output, input_file, stats_file)

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
                
                    try:
                        output = subprocess.check_output(command.format(os.path.abspath(input_file)).split())
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

                    except subprocess.CalledProcessError:
                        log.warn("{} failed running on {}".format(command, input_file))
                        raise SystemExit
            stats.write("\n\n")