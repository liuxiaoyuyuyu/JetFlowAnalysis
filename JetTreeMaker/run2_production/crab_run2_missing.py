"""CRAB configuration for the missing Run-2 UL QCD AK8 trees.

Select exactly one entry from missing_samples.json with RUN2_SAMPLE, for
example:

  RUN2_SAMPLE=2018_qcd_600to800 crab submit -c crab_run2_missing.py

Keeping the sample selection outside this file makes every CRAB task explicit
and prevents importing the configuration from submitting all five samples.
"""

import json
import os

from CRABClient.UserUtilities import config as CrabConfig


HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "missing_samples.json")) as handle:
    SAMPLES = {sample["key"]: sample for sample in json.load(handle)["samples"]}

sample_key = os.environ.get("RUN2_SAMPLE")
if sample_key not in SAMPLES:
    choices = ", ".join(sorted(SAMPLES))
    raise RuntimeError(
        "Set RUN2_SAMPLE to exactly one of: {}".format(choices)
    )

sample = SAMPLES[sample_key]
era_label = {
    "2016APV": "2016APV",
    "2016": "2016",
    "2017": "2017",
    "2018": "2018",
}[sample["era"]]

request_name = "Run2UL_{}_AK8Trees_revisit_v1".format(sample_key)
output_tag = (
    "JetTrees_PFJetsAK8_Pythia8_Pt{}_Run{}_MiniAODv2_20261111_revisit_v1"
    .format(sample["qhat"], era_label)
)

config = CrabConfig()

config.section_("General")
config.General.requestName = request_name
config.General.workArea = os.path.join(HERE, "crab_projects_run2_missing")
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = os.path.join(HERE, "run2_mc_cfg.py")
config.JobType.outputFiles = ["output_UL_ak8.root"]
config.JobType.pyCfgParams = [
    "globalTag={}".format(sample["global_tag"]),
    "outputFile=output_UL_ak8.root",
]
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 3000
config.JobType.maxJobRuntimeMin = 2750

config.section_("Data")
config.Data.inputDataset = sample["dataset"]
config.Data.inputDBS = "global"
# One MiniAOD file per job reproduces the granularity of the existing trees and
# keeps ROOT outputs comfortably below the CRAB transfer-size limit.
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.outputDatasetTag = output_tag
config.Data.outLFNDirBase = (
    "/store/group/phys_heavyions/xiaoyul/Run2_root_out/revisit/trees"
)

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
