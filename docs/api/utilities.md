# Utilities API

Helper functions and utility classes.

## System Information

### check_dependencies()

Check PRISM dependencies.

```python
import prism

deps = prism.check_dependencies()

for tool, available in deps.items():
    status = "✓" if available else "✗"
    print(f"{status} {tool}")
```

**Returns**: dict of dependency status

### list_forcefields()

List available GROMACS force fields.

```python
import prism

ffs = prism.list_forcefields()
for ff in ffs:
    print(ff)
```

### get_version()

Get PRISM version.

```python
import prism

version = prism.get_version()
print(f"PRISM version: {version}")
```

---

## Configuration

### ConfigurationManager

Manage configuration files.

```python
from prism.utils.config import ConfigurationManager

config = ConfigurationManager("config.yaml")
config.validate()
config.show()
```

---

## Environment

### GromacsEnvironment

GROMACS environment handling.

```python
from prism.utils.environment import GromacsEnvironment

env = GromacsEnvironment()
gmx_version = env.get_version()
gmx_path = env.get_gmx_command()
```

---

## Examples

```python
import prism

# Check setup
deps = prism.check_dependencies()
if not deps['gromacs']:
    print("ERROR: GROMACS not found")

# List force fields
ffs = prism.list_forcefields()
print(f"Available force fields: {len(ffs)}")

# Version info
print(f"PRISM version: {prism.get_version()}")
```

---

## Related

- [User Guide](../user-guide/index.md)
- [Installation](../getting-started/installation.md)
