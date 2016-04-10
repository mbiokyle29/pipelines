PATH_ROOT = "/sra/sra-instant/reads/ByExp/sra/SRX/SRX{}/{}/"

def read_srx_list(file):

    srxs = set()
    with open(file, "r") as fh:
        next(fh)
        for line in fh:
            if line.startswith("SRX"):
                srxs.add(line.rstrip())

    return list(srxs)

def read_srr_file(srr_file):

    srr_data = {}
    with open(srr_file, "r") as fh:
        
        for line in fh:
            line = line.rstrip()
            srx, srrs = line.split("\t")
            srrs = srrs.split(", ")
            assert srx.startswith("SRX")
            assert srx not in srr_data
            for srr in srrs:
                assert srr.startswith("SRR")
                assert "," not in srr
                assert " " not in srr

            srr_data[srx] = srrs

    return srr_data

def fetch_srr_ids(srx, conn):

    assert srx.startswith("SRX")
    first_chunk = srx[3:6]
    path = PATH_ROOT.format(first_chunk, srx)
    res = conn.nlst(path)
    return map(lambda x: x.split("/")[-1], res)

def build_srr_list(srxs, conn):

    srr_data = {}

    for srx in srx_ids:
        srrs = fetch_srr_ids(srx, conn)
        srr_data[srx] = srrs

    return srr_data
