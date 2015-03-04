#

from pymatgen.io.abinitio.flows import phonon_flow

from dftscripts.tasks.faketasks import get_qpts

__all__ = ['auto_phonon_flow']

def auto_phonon_flow(workdir, manager, inp, ngqpt = (4,4,4), with_nscf = False):
    qpts, mem = get_qpts(inp.deepcopy(), ngqpt, workdir = workdir)
    ph_inputs = []
    for qpt in qpts:
        ph_input = inp.deepcopy()
        ph_input.set_vars(dict(qpt = qpt))
        ph_inputs.append(ph_input)
    return phonon_flow(workdir, manager, inp, ph_inputs, with_nscf = with_nscf)
#
#class PhononFlow(AbinitFlow):
#    def __init__(self, workdir, manager, inp, pickle_protocol = -1, ngqpt = (4,4,4), do_nscf = True):
#        super(PhononFlow, self).__init__(workdir, manager, pickle_protocol)
#        self.set_workdir(workdir)
#        self.creation_date = time.asctime()
#        self.manager = manager.deepcopy()
#        self.inp = inp
#        self._ngqpt = ngqpt
#        self.do_nscf = do_nscf
#        self.nscf_tolwfr = 1.e-15
#        self.ph_tolvrs = 1.e-12
#
#    @property
#    def ngqpt(self):
#        return self._ngqpt
#
##    def nscf_input(self, qpt):
##        inp = self.inp.deepcopy()
##        inp.set_vars(
##            kptopt = 3,
##            nqpt = 1,
##            qpt = qpt,
##            iscf = -2,
##            tolvrs = 0,
##            tolwfr = self.nscf_tolwfr,
##        )
##        return inp
##
##    def ph_input(self, qpt, idir, ipert):
##        rfdir = 3 * [0]
##        rfdir[idir - 1] = 1
##        rfatpol = [ipert, ipert]
##        inp = self.inp.deepcopy()
##        inp.set_vars(
##            kptopt = 3,
##            nqpt = 1,
##            qpt = qpt,
##            rfphon = 1,
##            rfdir = rfdir,
##            rfatpol = rfatpol,
##            tolwfr = 0,
##            tolvrs = self.ph_tolvrs,
##        )
##        return inp
##
##    def merge_ddb_files(self):
##        ddb_files = filter(None, [workflow.outdir.has_abiext("DDB") for workflow in self])
##        out_ddb = self.outdir.path_in("out_DDB")
##        desc = "DDB file merged by %s on %s" % (self.__class__.__name__, time.asctime())
##        mrgddb = wrappers.Mrgddb(verbose=1)
##        mrgddb.set_mpi_runner("")
##        mrgddb.merge(ddb_files, out_ddb = out_ddb, description = desc, cwd=self.outdir.path)
##        return out_ddb
##        
##    def anaddb(self, qpath):
##        ddb_out = self.outdir.path_in("out_DDB")
##        if not os.path.exists(ddb_out):
##            ddb_out = self.merge_ddb_files()
##            
##        anaddb = wrappers.Anaddb(verbose=1)
##        anaddb.set_mpi_runner("mpirun")
##        anaddb.diagonalize_qpath(ddb_out, self.ngqpt, qpath, cwd=self.outdir.path)
#
#    def build_tasks(self):
#        scf_input = self.inp.deepcopy()
#        qpts, mem = get_qpts(self.inp.deepcopy(), self.ngqpt, workdir = self.workdir)
#        works_ph = []
#        
#        for qpt in qpts:
#            if self.do_nscf:
#                nscf_input = self.nscf_input(qpt)
#                nscf_task = work_nscf.register(nscf_input, deps = {scf_task: ["WFK", "DEN"]}, task_class=NscfTask)
#            irred_perts, mem = get_irred_perts(qpt, self.inp.deepcopy(), manager = self.manager, workdir = self.workdir)
#            work_ph = PhononWorkflow()
#            print('qpt: ' + str(qpt) + ' ntasks: ' + str(len(irred_perts)))
#            for irred_pert in irred_perts:
#                idir = irred_pert["idir"]
#                ipert = irred_pert["ipert"]
#                ph_input = self.ph_input(qpt, idir, ipert)
#                if nscf_task is not None:
#                    work_ph.register(ph_input, deps={nscf_task: "WFQ", scf_task: "WFK"}, task_class=PhononTask)
#                else:
#                    work_ph.register(ph_input, deps={scf_task: "WFK"}, task_class=PhononTask)
#            works_ph.append(work_ph)
#        if work_nscf is not None:
#            self.register_work(work_nscf)
#
#        for work_ph in works_ph:
#            self.register_work(work_ph)
#        