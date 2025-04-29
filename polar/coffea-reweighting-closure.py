import uproot
import matplotlib.pyplot as plt

rootfiles = ['coffea-reweighting-test.root', 'reweight_for_ran_jhugen_0p3m_2e2mu.root']
trees = list(map(lambda x: uproot.open(x)['tree'].arrays(), rootfiles))
fields = trees[0].fields
print(fields)

for field in fields:
    plt.hist(trees[0][field], bins=100, density=True, histtype='step', label='Coffea')
    plt.hist(trees[1][field], bins=100, density=True, histtype='step', label='C++')
    plt.legend()
    plt.xlabel(field)
    plt.ylabel('Density')
    plt.grid()
    plt.tight_layout()
    plt.savefig('coffea-reweighting-closure-' + field + '.pdf')
    plt.clf()
