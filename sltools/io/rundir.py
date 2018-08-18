import os
import shutil
from datetime import datetime

class Rundir:
    def __init__(self, template):
        now = datetime.now().isoformat().replace(':', '')
        self.path = '%s-%s' % (template, now)
        os.mkdir(self.path)

    def _mirror(self, src, islink=False):
        name = os.path.basename(src)
        dst = '%s/%s' % (self.path, name)
        src = os.path.abspath(src)

        if islink:
            os.symlink(src, dst)
        else:
            shutil.copy2(src, dst)
        return name

    def link(self, src):
        return self._mirror(src, islink=True)

    def copy(self, src):
        return self._mirror(src)
