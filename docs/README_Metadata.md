As images are read via numpy files, voxelsize and other spatial information is missing.
To remedy this, metadata is read automatically from a .txt file with the same name.

This readme explaines the format we expect inside the txt file

# Comments

Comments are per line and start with #

# Key/Value Pairs

Each line might have one key/value pair. The key and value are seperated by either whitespace, a colon or an equal sign. The value part might contain multiple values. Different values must be seperated by either whitespace, comma or semicolon. They can also start and end with braces or brackets.

# Keys

| Keyword     | # values | Description                                  | Anchor | Default |
|-------------|----------|----------------------------------------------|--------|---------|
| bbox        | 2 or 6   | The minimum of all points, then the maximum. | corner | ---     |
| corner-bbox | 2 or 6   | The minimum of all points, then the maximum. | corner | ---     |
| center-bbox | 2 or 6   | The minimum of all points, then the maximum. | center | ---     |
| origin      | 1 or 3   | The minimum of all points.                   | center | 0.0     |
| corner      | 1 or 3   | The minimum of all points.                   | corner | ---     |
| center      | 1 or 3   | The minimum of all points.                   | center | 0.0     |
| voxelsize   | 1 or 3   | Voxelsize                                    | -----  | 1.0     |

# Priority

As one can calculate the bbox from the origin and voxelsize and vice versa, we have the following priority for any metadata file (from highest to lowest):

- Usage of origin/corner/center and voxelsize if both are defined.
- Usage of bbox if defined.
- Usage of origin/corner/center and voxelsize with default values if either (or even both) are not defined.