"""
Tests for convgeno.slurm.config

Run with:  pytest tests/test_slurm_config.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from convgeno.slurm.config import (
    MultinodeConfig,
    PhenotypeTreeConfig,
    PipelineConfig,
    SlurmConfig,
    SpeciesTreeConfig,
    UltrametricConfig,
    normalize_optional_account,
)


class TestUltrametricConfig:
    def test_calibration_cli_arg_and_has_calibration(self):
        cfg = UltrametricConfig(
            species_a="Homo_sapiens", species_b="Felis_catus", divergence_my=94.0
        )
        assert cfg.has_calibration() is True
        assert (
            cfg.calibration_cli_arg()
            == "Homo_sapiens_Felis_catus:Homo_sapiens,Felis_catus:94"
        )

    def test_empty_has_no_calibration(self):
        cfg = UltrametricConfig()
        assert cfg.has_calibration() is False
        assert cfg.calibration_cli_arg() is None

    def test_from_dict_ignores_unknown_keys(self):
        cfg = UltrametricConfig.from_dict(
            {"species_a": "a", "species_b": "b", "divergence_my": 10, "future": "x"}
        )
        assert cfg.species_a == "a" and cfg.divergence_my == 10

    def test_pipeline_config_round_trips_calibration(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            ultrametric=UltrametricConfig(
                species_a="Homo_sapiens", species_b="Felis_catus", divergence_my=94.0
            ),
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        assert "ultrametric:" in path.read_text(encoding="utf-8")
        loaded = PipelineConfig.load(path)
        assert loaded.ultrametric is not None
        assert loaded.ultrametric.species_a == "Homo_sapiens"
        assert loaded.ultrametric.divergence_my == 94.0

    def test_pipeline_config_omits_empty_calibration(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            ultrametric=UltrametricConfig(),  # skipped calibration
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        assert "ultrametric:" not in path.read_text(encoding="utf-8")
        assert PipelineConfig.load(path).ultrametric is None


class TestSpeciesTreeConfig:
    def test_has_tree(self):
        assert SpeciesTreeConfig(path="/t.nwk").has_tree() is True
        assert SpeciesTreeConfig().has_tree() is False

    def test_from_dict_ignores_unknown_keys(self):
        cfg = SpeciesTreeConfig.from_dict(
            {"path": "/t.nwk", "is_ultrametric": True, "num_sites": 500, "next": "x"}
        )
        assert cfg.path == "/t.nwk"
        assert cfg.is_ultrametric is True
        assert cfg.num_sites == 500

    def test_to_dict_omits_none_num_sites(self):
        d = SpeciesTreeConfig(path="/t.nwk").to_dict()
        assert d["path"] == "/t.nwk"
        assert d["is_ultrametric"] is False
        assert "num_sites" not in d  # None is omitted

    def test_pipeline_config_round_trips_species_tree(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            species_tree=SpeciesTreeConfig(
                path="/data/my_tree.nwk", is_ultrametric=False, num_sites=1234
            ),
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        assert "species_tree:" in path.read_text(encoding="utf-8")
        loaded = PipelineConfig.load(path)
        assert loaded.species_tree is not None
        assert loaded.species_tree.path == "/data/my_tree.nwk"
        assert loaded.species_tree.is_ultrametric is False
        assert loaded.species_tree.num_sites == 1234


class TestPhenotypeTreeConfig:
    def test_has_table(self):
        assert PhenotypeTreeConfig(table="/p.tsv").has_table() is True
        assert PhenotypeTreeConfig().has_table() is False

    def test_defaults(self):
        cfg = PhenotypeTreeConfig(table="/p.tsv")
        assert cfg.id_col == "species"
        assert cfg.pheno_col == "phenotype"
        assert cfg.model == "ER"

    def test_from_dict_ignores_unknown_keys(self):
        cfg = PhenotypeTreeConfig.from_dict(
            {
                "table": "/p.tsv",
                "id_col": "sp",
                "pheno_col": "diet",
                "model": "SYM",
                "next": "x",
            }
        )
        assert cfg.table == "/p.tsv"
        assert cfg.id_col == "sp"
        assert cfg.pheno_col == "diet"
        assert cfg.model == "SYM"

    def test_pipeline_config_round_trips_phenotype_tree(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            phenotype_tree=PhenotypeTreeConfig(
                table="/data/phenotypes.tsv",
                id_col="species",
                pheno_col="diet",
                model="SYM",
            ),
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        assert "phenotype_tree:" in path.read_text(encoding="utf-8")
        loaded = PipelineConfig.load(path)
        assert loaded.phenotype_tree is not None
        assert loaded.phenotype_tree.table == "/data/phenotypes.tsv"
        assert loaded.phenotype_tree.pheno_col == "diet"
        assert loaded.phenotype_tree.model == "SYM"

    def test_pipeline_config_omits_empty_phenotype_tree(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            phenotype_tree=PhenotypeTreeConfig(),  # no table
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        assert "phenotype_tree:" not in path.read_text(encoding="utf-8")
        assert PipelineConfig.load(path).phenotype_tree is None

    def test_pipeline_config_round_trips_ultrametric_tree(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            species_tree=SpeciesTreeConfig(path="/data/dated.nwk", is_ultrametric=True),
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.species_tree is not None
        assert loaded.species_tree.is_ultrametric is True
        assert loaded.species_tree.num_sites is None

    def test_pipeline_config_omits_when_no_tree(self, tmp_path: Path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            species_tree=SpeciesTreeConfig(),  # no path -> not written
        )
        path = tmp_path / "pipeline_config.yaml"
        cfg.save(path)
        assert "species_tree:" not in path.read_text(encoding="utf-8")
        assert PipelineConfig.load(path).species_tree is None


class TestSlurmConfigFromDict:
    def test_valid_dict(self):
        cfg = SlurmConfig.from_dict(
            {"partition": "hawkcpu", "cpus_per_task": 32, "time_limit": "24:00:00"}
        )
        assert cfg.partition == "hawkcpu"
        assert cfg.cpus_per_task == 32
        assert cfg.time_limit == "24:00:00"
        assert cfg.nodes == 1
        assert cfg.ntasks == 1
        assert cfg.mem is None
        assert cfg.mem_per_cpu is None
        assert cfg.account is None
        assert cfg.mail_user is None
        assert cfg.mail_type == "END,FAIL"
        assert cfg.extra_sbatch_args == []

    def test_missing_partition_raises(self):
        with pytest.raises(ValueError, match="partition"):
            SlurmConfig.from_dict({"cpus_per_task": 16})

    def test_ignores_unknown_keys(self):
        cfg = SlurmConfig.from_dict(
            {"partition": "rapids", "some_future_field": "value123"}
        )
        assert cfg.partition == "rapids"


class TestToSbatchLines:
    def test_basic(self):
        cfg = SlurmConfig(partition="hawkcpu")
        lines = cfg.to_sbatch_lines()
        assert isinstance(lines, list)
        assert "#SBATCH --partition=hawkcpu" in lines
        assert "#SBATCH --time=72:00:00" in lines
        assert "#SBATCH --cpus-per-task=16" in lines
        assert not any("--mem=" in line for line in lines)
        assert not any("--mem-per-cpu" in line for line in lines)
        assert "#SBATCH --nodes=1" in lines
        assert not any("--account" in line for line in lines)
        assert not any("--mail-user" in line for line in lines)

    def test_with_extra_args(self):
        cfg = SlurmConfig(
            partition="rapids",
            extra_sbatch_args=["--qos=high", "--reservation=myres"],
        )
        lines = cfg.to_sbatch_lines()
        assert "#SBATCH --qos=high" in lines
        assert "#SBATCH --reservation=myres" in lines

    def test_with_optional_fields(self):
        cfg = SlurmConfig(
            partition="hawkcpu",
            account="wym219",
            mail_user="abc@lehigh.edu",
        )
        lines = cfg.to_sbatch_lines()
        assert "#SBATCH --account=wym219" in lines
        assert "#SBATCH --mail-user=abc@lehigh.edu" in lines

    def test_mem_default_is_unset(self):
        cfg = SlurmConfig(partition="hawkcpu")

        assert not any("--mem=" in line for line in cfg.to_sbatch_lines())

    def test_mem_custom_value(self):
        cfg = SlurmConfig(partition="hawkcpu", mem="350400M")

        assert "#SBATCH --mem=350400M" in cfg.to_sbatch_lines()

    def test_mem_per_cpu_not_emitted_when_mem_is_set(self):
        cfg = SlurmConfig(partition="hawkcpu", mem="350400M")

        lines = cfg.to_sbatch_lines()

        assert "#SBATCH --mem=350400M" in lines
        assert not any("--mem-per-cpu" in line for line in lines)

    def test_mem_per_cpu_emitted_when_mem_is_unset(self):
        cfg = SlurmConfig(partition="hawkcpu", mem_per_cpu="4G")

        assert "#SBATCH --mem-per-cpu=4G" in cfg.to_sbatch_lines()

    def test_rejects_mem_and_mem_per_cpu_together(self):
        with pytest.raises(ValueError, match="mem.*mem_per_cpu"):
            SlurmConfig(partition="hawkcpu", mem="350400M", mem_per_cpu="4G")


class TestPipelineConfigRoundtrip:
    def test_save_and_load(self, tmp_path: Path):
        original = PipelineConfig(
            project_dir="/share/ceph/wym219group/shared/projects/GeneFamilyConverge",
            conda_env="convgeno",
            slurm=SlurmConfig(partition="hawkcpu", cpus_per_task=32),
        )
        config_path = tmp_path / "config.yaml"
        original.save(config_path)
        loaded = PipelineConfig.load(config_path)

        assert loaded.project_dir == original.project_dir
        assert loaded.conda_env == original.conda_env
        assert loaded.slurm.partition == "hawkcpu"
        assert loaded.slurm.cpus_per_task == 32
        assert loaded.slurm.time_limit == "72:00:00"
        assert loaded.slurm.nodes == 1

    def test_from_dict_backward_compat_mem_per_cpu(self):
        cfg = SlurmConfig.from_dict(
            {"partition": "hawkcpu", "mem_per_cpu": "4G"}
        )

        assert cfg.mem is None
        assert cfg.mem_per_cpu == "4G"

    def test_from_dict_with_mem_key(self):
        cfg = SlurmConfig.from_dict({"partition": "hawkcpu", "mem": "350400M"})

        assert cfg.mem == "350400M"
        assert cfg.mem_per_cpu is None

    def test_from_dict_rejects_mem_and_mem_per_cpu(self):
        with pytest.raises(ValueError, match="mem.*mem_per_cpu"):
            SlurmConfig.from_dict(
                {"partition": "hawkcpu", "mem": "350400M", "mem_per_cpu": "4G"}
            )

    def test_load_missing_file(self):
        with pytest.raises(FileNotFoundError):
            PipelineConfig.load("/nonexistent/path/config.yaml")


class TestToDictExcludesNone:
    def test_none_fields_omitted(self):
        cfg = SlurmConfig(partition="hawkcpu")
        d = cfg.to_dict()
        assert "account" not in d
        assert "mail_user" not in d
        assert "partition" in d


class TestOptionalAccountHandling:
    @pytest.mark.parametrize("account", [None, "", " ", "null", "None", "NULL"])
    def test_optional_account_values_are_omitted(self, account):
        cfg = SlurmConfig(partition="hawkcpu", account=account)
        lines = cfg.to_sbatch_lines()
        assert not any("--account" in line for line in lines)

    def test_valid_account_is_included(self):
        cfg = SlurmConfig(partition="hawkcpu", account="valid_alloc")
        lines = cfg.to_sbatch_lines()
        assert "#SBATCH --account=valid_alloc" in lines

    @pytest.mark.parametrize("account", [None, "", "null", "None"])
    def test_from_dict_normalizes_null_variants(self, account):
        data = {"partition": "hawkcpu"}
        if account is not None:
            data["account"] = account
        cfg = SlurmConfig.from_dict(data)
        assert cfg.account is None

    def test_from_dict_preserves_valid_account(self):
        cfg = SlurmConfig.from_dict({"partition": "hawkcpu", "account": "wym219"})
        assert cfg.account == "wym219"

    def test_to_dict_omits_normalized_null_account(self):
        cfg = SlurmConfig(partition="hawkcpu", account="null")
        lines = cfg.to_sbatch_lines()
        assert not any("--account" in line for line in lines)

    def test_save_load_roundtrip_null_account(self, tmp_path):
        config = PipelineConfig(
            project_dir="/project",
            slurm=SlurmConfig(partition="hawkcpu", account=None),
        )
        path = tmp_path / "config.yaml"
        config.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.slurm.account is None
        assert not any("--account" in line for line in loaded.slurm.to_sbatch_lines())


class TestScratchConfig:
    def test_scratch_dir_in_to_dict(self):
        cfg = SlurmConfig(partition="hawkcpu", scratch_dir="/scratch/user")

        assert cfg.to_dict()["scratch_dir"] == "/scratch/user"

    def test_scratch_dir_default_none(self):
        cfg = SlurmConfig(partition="hawkcpu")

        assert "scratch_dir" in cfg.to_dict()
        assert cfg.to_dict()["scratch_dir"] is None

    def test_scratch_dir_from_dict_missing(self):
        cfg = SlurmConfig.from_dict({"partition": "hawkcpu"})

        assert cfg.scratch_dir is None

    def test_scratch_dir_not_in_sbatch_lines(self):
        cfg = SlurmConfig(partition="hawkcpu", scratch_dir="/scratch/user")

        assert not any("scratch" in line for line in cfg.to_sbatch_lines())

    def test_is_ephemeral_scratch_roundtrip(self):
        original = SlurmConfig(
            partition="hawkcpu",
            scratch_dir="/local/scratch",
            is_ephemeral_scratch=True,
        )

        restored = SlurmConfig.from_dict(original.to_dict())

        assert restored.scratch_dir == "/local/scratch"
        assert restored.is_ephemeral_scratch is True


class TestMultinodeConfig:
    def test_defaults_all_none_and_empty_dict(self):
        mn = MultinodeConfig()
        assert mn.throughput_const is None
        assert mn.to_dict() == {}  # None fields omitted

    def test_to_dict_omits_none(self):
        mn = MultinodeConfig(waves=6, throughput_const=1.23e12)
        assert mn.to_dict() == {"waves": 6, "throughput_const": 1.23e12}

    def test_from_dict_ignores_unknown_keys(self):
        mn = MultinodeConfig.from_dict(
            {"waves": 3, "search_cpus": 20, "future_knob": "ignored"}
        )
        assert mn.waves == 3
        assert mn.search_cpus == 20

    def test_pipeline_config_roundtrip_with_multinode(self, tmp_path):
        cfg = PipelineConfig(
            project_dir="/p",
            slurm=SlurmConfig(partition="hawkcpu"),
            multinode=MultinodeConfig(
                array_throttle=12, throughput_const=5e11, search_mem="24000M"
            ),
        )
        path = tmp_path / "config.yaml"
        cfg.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.multinode is not None
        assert loaded.multinode.array_throttle == 12
        assert loaded.multinode.throughput_const == 5e11
        assert loaded.multinode.search_mem == "24000M"

    def test_pipeline_config_roundtrip_without_multinode(self, tmp_path):
        cfg = PipelineConfig(
            project_dir="/p", slurm=SlurmConfig(partition="hawkcpu")
        )
        path = tmp_path / "config.yaml"
        cfg.save(path)
        loaded = PipelineConfig.load(path)
        assert loaded.multinode is None
