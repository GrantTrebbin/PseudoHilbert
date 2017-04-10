"""
Generate Pseudo-Hilbert Curves for arbitrary rectangular areas.

Implementation of the paper
    "A Pseudo-Hilbert Scan for Arbitrarily-Sized Arrays" by
    Jian ZHANG, Sei-ichiro KAMATO, Yoshifumi UESHIGE

An arbitrary rectangular region is subdivided into blocks.  There are n x n
blocks where n is a power of two.  Blocks contain cells.  A hilbert curve is
used to traverse the blocks in a particular order.  The cells in each block
are further traverse in a manner so that entry and exit points of adjacent
cells line up.

Basic sub-types of the Hilbert curve

    Type  1            Type  2            Type  3            Type  4
(0,1)     (1,1)    (0,1)     (1,1)    (0,1)     (1,1)    (0,1)     (1,1)
  *---------*        <---------*        ^         v        *---------<
  |         |                  |        |         |        |
  |         |                  |        |         |        |
  |         |                  |        |         |        |
  ^         v        >---------*        *---------*        *--------->
(0,0)     (1,0)    (0,0)     (1,0)    (0,0)     (1,0)    (0,0)     (1,0)


Directions are defined as follows

^ UP     +y
v DOWN   -y
> RIGHT  +x
< LEFT   -x

Coordinate system
 +y
 |
 |
 o----+x

(x,y)
unless specified otherwise, coordinates follow a (width, height) pattern
"""
import math
import operator
from enum import Enum, auto
from itertools import accumulate


class Direction(Enum):
    """Track directions."""

    UP = auto()  # +y
    DOWN = auto()  # -y
    LEFT = auto()  # -x
    RIGHT = auto()  # +x


class Parity(Enum):
    """Tracks the oddness or evenness of an item."""

    EVEN = auto()
    ODD = auto()


class Block:
    """A Block is a region of cells to be scanned in a specific manner.

    Args:
        hilbert_type (int): Can be 1, 2, 3, 4. Defines the shape of sub-curves
        address_x (list): A list of ones and zeros that defines the block's
                          x location in a hilbert curve
        address_y (list): A list of ones and zeros that defines the block's
                          y location in a hilbert curve

    Attributes:
        address_x (list): Addresses are a binary list of ones and zeros.
                          Most significant bit is in the zero position
        address_y (list): Addresses are a binary list of ones and zeros.
                          Most significant bit is in the zero position

        x_index (int): Decimal representation of the binary addresses
        y_index (int): Decimal representation of the binary addresses

        travel_direction_to_enter (Direction): Directions into a block
        travel_direction_to_leave (Direction): Directions out of a block

        x_size (int): The width of a block
        y_size (int): The height of a block

        x_pos (int): x position of the bottom left corner
        y_pos (int): y position of the bottom left corner

        shape (tuple): a two element tuple that holes the parity (odd or even)
                       of the width followed by the height of the block

        scan_type (int): Used to access information about where to start a
                         bidirectional raster scan of the block and which way
                         to go first
    """

    # Scan order instructions are a tuple three elements long
    #
    # 1 <= scan_type <= 8
    #
    # Zeroth element: if 0 start on the left side of the block
    #                if 1 start on the right side of the block
    # First element: if 0 start on the bottom side of the block
    #                if 1 start on the top side of the block
    # Second element: if 1 scan in the x direction and then the y direction
    #                 if 0 scan in the y direction and then the x direction

    # As the paper uses one-based indexing for the scanning types, 'None' has
    # been added in the zero position
    scan_instructions = [None,
                         (0, 0, 0),
                         (0, 0, 1),
                         (1, 1, 0),
                         (1, 1, 1),
                         (1, 0, 0),
                         (0, 1, 1),
                         (0, 1, 0),
                         (1, 0, 1)]

    # An improvement not included in the paper but added later by the author
    # changes the scan pattern for blocks with even-even parity. With both side
    # lengths greater than or equal to four.  The scan in broken down into
    # sub-blocks.  For example a block with a scan type of 7 can be replaced
    # with 4 blocks of type 6,7,7,8 with relative positions
    # (0, 1), (0, 0), (1, 0), (1, 1).  Each tuple below contains 4 tuples.
    # Each of these hold (Type, relative_x_position, relative_y_position)
    even_even_optimisation_path = [None,
                                   ((2, 0, 0), (1, 0, 1), (1, 1, 1), (4, 1, 0)),
                                   ((1, 0, 0), (2, 1, 0), (2, 1, 1), (3, 0, 1)),
                                   ((4, 1, 1), (3, 1, 0), (3, 0, 0), (2, 0, 1)),
                                   ((3, 1, 1), (4, 0, 1), (4, 0, 0), (1, 1, 0)),
                                   ((8, 1, 0), (5, 1, 1), (5, 0, 1), (6, 0, 0)),
                                   ((7, 0, 1), (6, 1, 1), (6, 1, 0), (5, 0, 0)),
                                   ((6, 0, 1), (7, 0, 0), (7, 1, 0), (8, 1, 1)),
                                   ((5, 1, 0), (8, 0, 0), (8, 0, 1), (7, 1, 1))]

    def __init__(self, hilbert_type, address_x, address_y):
        """Initialize a block with its type and location in the Hilbert curve.

        Args
            hilbert_type (int): Can be 1, 2, 3, 4. Defines the
                                shape of sub-curves
            address_x (list): A list of ones and zeros that defines the block's
                              x location in a hilbert curve
            address_y (list): A list of ones and zeros that defines the block's
                              x location in a hilbert curve
        """
        self.hilbert_type = hilbert_type

        # Addresses are a binary list of ones and zeros
        # Most significant bit is in the zero position
        self.address_x = list(address_x)
        self.address_y = list(address_y)

        # Decimal representation of the binary addresses
        # The calculation of these are purposefully delayed and are set by
        # calling the calculate_decimal_indices method only after the block is
        # in it's final position
        self.x_index = None
        self.y_index = None

        self.travel_direction_to_enter = None
        self.travel_direction_to_leave = None

        self.x_size = None
        self.y_size = None

        self.x_pos = None
        self.y_pos = None

        self.shape = None

        # Holds information about where to start a bidirectional raster scan of
        # the block and which way to go first
        self.scan_type = None

    def calculate_decimal_indices(self):
        """Turn the binary hilbert indices into decimal integers."""
        x_bits = self.address_x
        x = 0
        for bit in x_bits:
            x = (x << 1) | bit
        self.x_index = x

        y_bits = self.address_y
        y = 0
        for bit in y_bits:
            y = (y << 1) | bit
        self.y_index = y

    def position_block(self, parent_block):
        """Position a Block within a parent Block.

           Set the binary address of this block to the binary address of this
           block appended to the address of a parent block.

        Args
            parent_block (block): A Block that this Block will be positioned in
        """
        self.address_x =\
            operator.add(parent_block.address_x, self.address_x)

        self.address_y =\
            operator.add(parent_block.address_y, self.address_y)

    def set_coordinates(self, x_coordinate, y_coordinate):
        """Set the location of the bottom left cell of the Block.

        Args
            x_coordinate (int): x position of the bottom left cell of the Block
            y_coordinate (int): y position of the bottom left cell of the Block
        """
        self.x_pos = x_coordinate
        self.y_pos = y_coordinate

    def set_size(self, x_size, y_size):
        """Set the width and height of the Block.

        Args
            x_size (int):
            y_size (int):
        """
        self.x_size = x_size
        self.y_size = y_size
        self.shape = (Parity.EVEN if x_size % 2 == 0 else Parity.ODD,
                      Parity.EVEN if
                      y_size % 2 == 0 else Parity.ODD)

    def copy(self):
        """Return a copy of the Block object.

        Returns
            new_block (Block): A new Block with attributes set
        """
        new_block = Block(self.hilbert_type, self.address_x, self.address_y)
        new_block.travel_direction_to_enter = self.travel_direction_to_enter
        new_block.travel_direction_to_leave = self.travel_direction_to_leave
        new_block.x_index = self.x_index
        new_block.y_index = self.y_index
        new_block.x_size = self.x_size
        new_block.y_size = self.y_size
        new_block.x_pos = self.x_pos
        new_block.y_pos = self.y_pos
        return new_block

    def scan(self):
        """Determine how to scan the Block.

        If both sides are divisible by 4, do a further division. If not, go
        straight to a bidirectional raster scan.

        Returns
            coordinates (list): A list of two element lists containing (x,y)
                                coordinates
        """
        coordinates = []
        if (self.x_size % 4 == 0 and
                self.y_size % 4 == 0):
            half_x_size = self.x_size // 2
            half_y_size = self.y_size // 2

            # Create the sub-blocks
            sub_blocks = [Block(None, [], []),
                          Block(None, [], []),
                          Block(None, [], []),
                          Block(None, [], [])]

            # Set the parameters for the sub-blocks
            for block_index in range(4):
                sub_blocks[block_index].set_size(half_x_size, half_y_size)

                sub_blocks[block_index].set_coordinates(
                    self.x_pos +
                    self.even_even_optimisation_path[
                        self.scan_type][block_index][1] * half_x_size,
                    self.y_pos +
                    self.even_even_optimisation_path[
                        self.scan_type][block_index][2] * half_y_size)

                sub_blocks[block_index].scan_type =\
                    self.even_even_optimisation_path[
                        self.scan_type][block_index][0]
                coordinates = coordinates +\
                    sub_blocks[block_index].bidirectional_raster_scan()
        else:
            coordinates = self.bidirectional_raster_scan()
        return coordinates

    def bidirectional_raster_scan(self):
        """Perform a bidirectional raster scan of a Block.

        Starts at a specified corner and heads in the x or y direction first.

        Returns
            coordinates (list): A list of two element lists containing (x,y)
                                coordinates
        """
        #  scan        scan       scan        scan
        #  Type 1      Type 2     Type 3      Type 4
        #  ---    ^    ------>    ---    o    ------o
        # |   |   |   |          |   |   |   |
        # |   |   |    ------    |   |   |    ------
        # |   |   |          |   |   |   |          |
        # o    ---    o------    v    ---    <------

        #  scan        scan       scan        scan
        #  Type 5      Type 6     Type 7      Type 8
        # ^    ---    o------    o    ---    <------
        # |   |   |          |   |   |   |          |
        # |   |   |    ------    |   |   |    ------
        # |   |   |   |          |   |   |   |
        #  ---    o    ------>    ---    v    ------o

        # Scan order instructions a tuple three elements long
        # 1 <= scan_type <= 8
        # Zeroth element: if 0 start on the bottom side of the block
        #                if 1 start on the top side of the block
        # First element: if 0 start on the left side of the block
        #                if 1 start on the right side of the block
        # Second element: if 1 scan in the x direction and then the y direction
        #                 if 0 scan in the y direction and then the x direction
        instructions = self.scan_instructions[self.scan_type]

        # Find the coordinates of the first cell of the block
        x = self.x_pos + instructions[0] * (self.x_size - 1)
        y = self.y_pos + instructions[1] * (self.y_size - 1)

        # Create two list "vectors" of which way to travel in each direction
        # Assumes scanning in the y direction first
        # primary equals [0, -1] or [0, 1] secondary equals [-1, 0] or [1, 0]
        primary_scan_direction = [0, 1 - 2 * instructions[1]]
        secondary_scan_directions = [1 - 2 * instructions[0], 0]

        # Size in the primary and secondary directions
        primary_size = self.y_size
        secondary_size = self.x_size

        # To scan in the x direction and then the y direction swap the
        # directions sizes and the scan directions
        if instructions[2] == 1:
            primary_scan_direction, secondary_scan_directions =\
                secondary_scan_directions, primary_scan_direction
            primary_size, secondary_size = secondary_size, primary_size

        # Go backwards one cell before starting the loop as the first
        # operation is a move forward.  This makes the first cell coordinate
        # correct and makes the loop cleaner
        x -= secondary_scan_directions[0]
        y -= secondary_scan_directions[1]

        # Initialize an empty list of coordinates for the cell
        coordinates = [None] * (self.x_size * self.y_size)

        counter = 0
        for secondary_counter in range(secondary_size):
            # Move in the secondary scan direction
            x += secondary_scan_directions[0]
            y += secondary_scan_directions[1]
            coordinates[counter] = [x, y]
            counter += 1

            for primary_counter in range(primary_size - 1):
                # Move in the primary scan direction
                x += primary_scan_direction[0]
                y += primary_scan_direction[1]
                coordinates[counter] = [x, y]
                counter += 1

            # Reverse the primary scan direction
            primary_scan_direction[0] *= -1
            primary_scan_direction[1] *= -1

        return coordinates


class PseudoHilbert:
    """Hold information about and generate pseudo Hilbert curves.

    Generate a pseudo Hilbert curve for an arbitrary rectangular region
    The curve will start in the lower left cell of the region and finish
    somewhere close to the bottom right cell visiting every cell

    Args
        width (int): Width of the arbitrary rectangle
        height (int): Height of the arbitrary rectangle

    Attributes
        width (int):
        height (int):
        coordinate_to_index (list): A 2D list of lists to look up the index from
                                    coordinates
        index_to_coordinate (list): A list of two element (x,y) list coordinates
                                    to lookup the coordinates from an index

    """

    # Lookup table for directions when given two consecutive Hilbert block types
    # Only works for blocks with EVEN EVEN parity
    even_even_block_directions = [None,
                                  [None,
                                   Direction.RIGHT,
                                   Direction.RIGHT,
                                   Direction.DOWN,
                                   Direction.DOWN],
                                  [None,
                                   Direction.UP,
                                   Direction.UP,
                                   Direction.LEFT,
                                   Direction.LEFT],
                                  [None,
                                   Direction.UP,
                                   Direction.UP,
                                   Direction.LEFT,
                                   Direction.LEFT],
                                  [None,
                                   Direction.RIGHT,
                                   Direction.RIGHT,
                                   Direction.DOWN,
                                   Direction.DOWN]]

    # When generating the hilbert curve an initial block is replaced by
    # 4 sub blocks.  This process is repeated until the required order is
    # reached.  This lookup table holds the subtypes and shape for each block.
    # 'None' is used because of the one based indexing used in the paper.
    template_table = [None,
                      [Block(2, [0], [0]), Block(1, [0], [1]),
                       Block(1, [1], [1]), Block(4, [1], [0])],
                      [Block(1, [0], [0]), Block(2, [1], [0]),
                       Block(2, [1], [1]), Block(3, [0], [1])],
                      [Block(4, [1], [1]), Block(3, [1], [0]),
                       Block(3, [0], [0]), Block(2, [0], [1])],
                      [Block(3, [1], [1]), Block(4, [0], [1]),
                       Block(4, [0], [0]), Block(1, [1], [0])]]

    # When at least one of the edges of the overall region is odd,
    # these are the scan types to use for even-even blocks.
    # For example, if you enter a block by going right and leave it by going up
    # the scan type is 7
    either_odd_scan_lookup = {(Direction.LEFT, Direction.LEFT): 5,
                              (Direction.LEFT, Direction.UP): 8,
                              (Direction.LEFT, Direction.DOWN): 5,
                              (Direction.RIGHT, Direction.RIGHT): 7,
                              (Direction.RIGHT, Direction.UP): 7,
                              (Direction.RIGHT, Direction.DOWN): 6,
                              (Direction.UP, Direction.LEFT): 5,
                              (Direction.UP, Direction.RIGHT): 8,
                              (Direction.UP, Direction.UP): 8,
                              (Direction.DOWN, Direction.LEFT): 6,
                              (Direction.DOWN, Direction.RIGHT): 7,
                              (Direction.DOWN, Direction.DOWN): 6,
                              (None, Direction.UP): 0,
                              (None, Direction.RIGHT): 0,
                              (Direction.RIGHT, None): 7,  # The last block
                              (Direction.DOWN, None): 7}  # The last block

    @staticmethod
    def division(length):
        """Divide a length into two sub lengths.

        This is used to determine the size of the blocks


        Args:
            length (int): the length of a block side to divide.
        Returns:
            (list): contains the first and second part of the divided length.
        """
        m = math.frexp(length)[1] - 2

        two_to_the_mth = pow(2, m)
        cutoff = 3 * two_to_the_mth
        second_part = two_to_the_mth

        if length > cutoff:
            second_part = two_to_the_mth * 2

        first_part = length - second_part
        return [first_part, second_part]

    def hilbert_type_to_direction(self, block_1, block_2):
        """Return the direction of travel when going from one block to another.

        Works for even-even blocks only.

        Args:
            block_1 (Block): The Block the path begins on
            block_2 (Block): The Block the path ends on

        Returns:
            (Direction):
        """
        # When travelling between two blocks the direction can be determined
        # by the hilbert type of each block.  This performs that lookup.
        return (self.even_even_block_directions[block_1.hilbert_type]
                )[block_2.hilbert_type]

    @staticmethod
    def set_scan_directions_even_even(block_list):
        """Set the scan direction of blocks in an even-even region.

        When the arbitrary rectangle is of type even-even the scan type of each
        block is equal to its hilbert type

        Args:
            block_list (list): A list of blocks that need their scan type set.
        """
        for block in block_list:  # type: Block
            block.scan_type = block.hilbert_type

    def set_scan_directions_either_odd(self, block_list):
        """Set the scan direction of blocks in non even-even regions.

        If the length of at least one edge of the arbitrary rectangle is odd,
        the block scan types need to bet set in a specific way.

        Args:
            block_list (list): A list of blocks that need their scan type set.
        """
        for block in block_list:  # type: Block

            # Set the scan type assuming that the shape is even-even
            block.scan_type =\
                self.either_odd_scan_lookup[(block.travel_direction_to_enter,
                                             block.travel_direction_to_leave)]

            # Fix the scan types for blocks that aren't even-even
            # These are predetermined values from the paper
            if block.shape == (Parity.ODD, Parity.EVEN):
                block.scan_type = 8
            if block.shape == (Parity.EVEN, Parity.ODD):
                block.scan_type = 7

        # The first block is a special case that has set scan types
        # depending on the shape
        if block_list[0].shape == (Parity.ODD, Parity.EVEN):
            block_list[0].scan_type = 1
        if block_list[0].shape == (Parity.EVEN, Parity.ODD):
            block_list[0].scan_type = 2
        if block_list[0].shape == (Parity.ODD, Parity.ODD):
            block_list[0].scan_type = 1

    def __init__(self, width, height):
        """Initialise and generate a Pseudo Hilbert Curve.

        Args:
            width (int): The width of the arbitrary rectangular region
            height (int): The height of the arbitrary rectangular region
        """
        self.width = width
        self.height = height
        self.coordinate_to_index = None
        self.index_to_coordinate = None

        # Minimum of floor(log2(width/2)) and floor(log2(height/2))
        # Determines the order of the parent Hilbert curve
        self.order = min(math.frexp(self.width)[1],
                         math.frexp(self.height)[1]) - 2

        x_divisions = [self.width]
        y_divisions = [self.height]

        # Calculate how to divide the arbitrary rectangle into blocks
        for order_count in range(self.order):
            for division_count in range(len(x_divisions)):
                length_to_divide = x_divisions.pop(0)
                divided_length = self.division(length_to_divide)
                x_divisions.extend(divided_length)

                length_to_divide = y_divisions.pop(0)
                divided_length = self.division(length_to_divide)
                y_divisions.extend(divided_length)

        # Calculate the coordinates for the lower left corner of the blocks
        cumulative_x_divisions = list(accumulate(x_divisions))
        cumulative_y_divisions = list(accumulate(y_divisions))
        cumulative_x_divisions.insert(0, 0)
        cumulative_x_divisions.pop()
        cumulative_y_divisions.insert(0, 0)
        cumulative_y_divisions.pop()

        # Generate the Hilbert curve in an iterative manner
        block_list_size = pow(2, 2*self.order)
        block_list = [Block] * block_list_size
        block_list[0] = Block(1, [], [])
        current_block_count = 1
        for order_count in range(self.order):

            # Replace blocks with sub-blocks
            for block_count in range(current_block_count-1, -1, -1):
                block = block_list[block_count]  # type: Block
                hilbert_type = block.hilbert_type

                # Create new block by copying them
                new_blocks = [block.copy() for block in
                              self.template_table[hilbert_type]]

                # Position new blocks within parent block
                for new_block in new_blocks:
                    new_block.position_block(block)

                # Add the blocks to the block list
                block_list[block_count * 4 + 0] = new_blocks[0]
                block_list[block_count * 4 + 1] = new_blocks[1]
                block_list[block_count * 4 + 2] = new_blocks[2]
                block_list[block_count * 4 + 3] = new_blocks[3]

            current_block_count *= 4

        # Set block travel directions
        for block_index in range(0, len(block_list) - 1):
            direction =\
                self.hilbert_type_to_direction(block_list[block_index],
                                               block_list[block_index + 1])

            first_block = block_list[block_index]
            second_block = block_list[block_index + 1]

            first_block.travel_direction_to_leave = direction
            second_block.travel_direction_to_enter = direction

        # Set the coordinates of each block
        for block in block_list:  # type: Block
            block.calculate_decimal_indices()
            block.set_size(x_divisions[block.x_index],
                           y_divisions[block.y_index])
            block.set_coordinates(cumulative_x_divisions[block.x_index],
                                  cumulative_y_divisions[block.y_index])

        # The shape of the overall arbitrary rectangle is the same as the first
        # block as all other row and column dimensions are even.
        first_block = block_list[0]  # type: Block
        overall_shape = first_block.shape
        if overall_shape == (Parity.EVEN, Parity.EVEN):
            self.set_scan_directions_even_even(block_list)
        else:
            self.set_scan_directions_either_odd(block_list)

        # Create the output lists.  i and j are unused.
        self.coordinate_to_index =\
            [[None for i in range(self.height)] for j in range(self.width)]
        self.index_to_coordinate =\
            [None for i in range(self.height * self.width)]

        # Scan each block and then fill the output lists
        counter = 0
        for block in block_list:  # type: Block
            scan_result = block.scan()
            for cell in scan_result:  # type: list
                self.coordinate_to_index[cell[0]][cell[1]] = counter
                self.index_to_coordinate[counter] = cell
                counter += 1

            # debugging code for blocks
            if False:
                print(block.hilbert_type,
                      block.address_x,
                      block.x_index,
                      block.address_y,
                      block.y_index,
                      block.x_pos,
                      block.y_pos,
                      block.x_size,
                      block.y_size,
                      block.travel_direction_to_enter,
                      block.travel_direction_to_leave,
                      block.shape,
                      block.scan_type)