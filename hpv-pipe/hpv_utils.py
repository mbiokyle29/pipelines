
def build_srx_list(file):

    srxs = set()
    with open(file, "r") as fh:
        for line in fh:
            if line.startswith("SRX"):
                srxs.add(line.rstrip())

    return list(srxs)

def fetch_srr_ids(srx, conn):

    assert srx.startswith("SRX")
    