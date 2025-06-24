"""
Default configuration settings for the project.
"""

from dataclasses import dataclass
from typing import Optional, List, Dict, Any

@dataclass
class ModelConfig:
    model_type: str = "transformer"  # or "diffusion"
    hidden_size: int = 768
    num_layers: int = 12
    num_heads: int = 12
    dropout: float = 0.1
    vocab_size: int = 32000

@dataclass
class TrainingConfig:
    batch_size: int = 32
    learning_rate: float = 1e-4
    num_epochs: int = 100
    warmup_steps: int = 1000
    gradient_accumulation_steps: int = 1
    max_grad_norm: float = 1.0
    weight_decay: float = 0.01
    device: str = "cuda"  # or "cpu"
    mixed_precision: bool = True

@dataclass
class DataConfig:
    train_data_path: str = "data/train"
    val_data_path: str = "data/val"
    test_data_path: str = "data/test"
    max_seq_length: int = 512
    num_workers: int = 4

@dataclass
class Config:
    model: ModelConfig = ModelConfig()
    training: TrainingConfig = TrainingConfig()
    data: DataConfig = DataConfig()
    seed: int = 42
    output_dir: str = "outputs"
    wandb_project: Optional[str] = "bioinformatics-dl"
    wandb_entity: Optional[str] = None
    checkpoint_dir: str = "checkpoints"
    log_dir: str = "logs" 