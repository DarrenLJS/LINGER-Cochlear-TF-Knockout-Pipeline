"""Make silent SGE-plugin submission errors visible. Usage: python patch_sge_logging.py [--revert]"""
import sys, shutil, importlib.util, pathlib
spec = importlib.util.find_spec("snakemake_executor_plugin_sge")
p = pathlib.Path(spec.origin)
bak = p.with_suffix(".py.orig")
if "--revert" in sys.argv:
    shutil.copy(bak, p); print("reverted", p); sys.exit()
src = p.read_text()
if "_sub_logged" in src:
    print("already patched"); sys.exit()
n = src.count("self._job_submission_executor.submit(")
assert n, "pattern not found - plugin version differs"
shutil.copy(p, bak)
src = src.replace("self._job_submission_executor.submit(", "self._sub_logged(")
helper = '''    def _sub_logged(self, fn, *args):
        """Submit to the thread pool and LOG any exception instead of swallowing it."""
        import traceback
        fut = self._job_submission_executor.submit(fn, *args)
        def _cb(f, _name=getattr(fn, "__name__", "submit")):
            exc = f.exception()
            if exc is not None:
                self.logger.error("SGE submission thread crashed in %s:\\n%s" % (
                    _name, "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))))
        fut.add_done_callback(_cb)
        return fut

'''
marker = "    def run_array_job(self"
assert marker in src
src = src.replace(marker, helper + marker, 1)
p.write_text(src)
print(f"patched {n} call(s) in {p}; backup at {bak}")
