import glob, collections
from concurrent.futures import ThreadPoolExecutor
files = sorted(glob.glob("/xdisk/clairemcwhite/proFAB/GO_sets/go_dataset_part*/*/random_positive_train_max1000.span")) + \
        sorted(glob.glob("/xdisk/clairemcwhite/proFAB/EC_sets/ec_dataset_part*/ecNo_*/random_positive_train_max1000.span"))
def rd(f):
    L = [l.split() for l in open(f) if l.strip()]
    return f, [x[0] for x in L], sum(len(x) > 1 for x in L)
with ThreadPoolExecutor(32) as ex: R = list(ex.map(rd, files))
tot = collections.Counter(); uniq = {"GO": set(), "EC": set()}; withspan = 0
for f, ids, ns in R:
    lib = "GO" if "GO_sets" in f else "EC"; tot[lib] += len(ids); uniq[lib] |= set(ids); withspan += ns
print({k: (tot[k], len(uniq[k])) for k in tot}, "all unique", len(uniq["GO"] | uniq["EC"]), "rows with explicit span coords", withspan, "files", len(files))
