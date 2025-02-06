from typing import Optional
from typing import Dict,List
class Section:
    """
    Represents a section with optional upper and lower identifiers.

    Attributes:
        upper (int | None): Identifier for the upper neighboring section.
        lower (int | None): Identifier for the lower neighboring section.
        length (float | None): Length of the section.
    """

    def __init__(self, upper: Optional[int] = None, lower: Optional[int] = None, length: Optional[float] = None):
        self.upper = upper
        self.lower = lower
        self.length = length

class SectionManager:
    """
    Manages a collection of Section objects and finds their neighbors.
    """

    def __init__(self, sections: List[Section]):
        self.sections = sections
        self.neighbours = self.find_neighbours()

    def find_neighbours(self) -> Dict[Section, List[Section]]:
        """
        Identifies neighboring sections based on their upper and lower identifiers.

        Returns:
            dict[Section, list[Section]]: A dictionary mapping each section to its neighboring sections.
        """
        # Create a mapping from identifiers to Section objects
        id_to_section = {}
        for section in self.sections:
            if section.upper is not None:
                id_to_section[section.upper] = section
            if section.lower is not None:
                id_to_section[section.lower] = section

        # Initialize the dictionary to hold neighbors
        neighbours = {section: [] for section in self.sections}

        # Populate the neighbors for each section
        for section in self.sections:
            if section.upper is not None and section.upper in id_to_section:
                neighbours[section].append(id_to_section[section.upper])
            if section.lower is not None and section.lower in id_to_section:
                neighbours[section].append(id_to_section[section.lower])

        return neighbours
    
    @staticmethod
    def generate_test_data(total: int = 1000) -> List[Section]:
        """
        Generates a list of Section objects for testing purposes.

        Args:
            total (int): The number of base sections to generate. Each base section will have
                         two additional sections with modified identifiers.

        Returns:
            list[Section]: A list containing the generated Section objects.
        """
        sections = []
        for idx in range(total):
            sections.append(Section(upper=None if idx == 0 else idx, lower=idx))
            sections.append(Section(upper=None if idx == 0 else idx + 2 * total, lower=idx + 2 * total + 1))
            sections.append(Section(upper=None if idx == 0 else idx + 4 * total, lower=idx + 4 * total + 1))
        return sections

# Generate test data
test_sections = SectionManager.generate_test_data(total=1000)

# Initialize the SectionManager with the generated test data
manager = SectionManager(test_sections)

# Access the generated sections
for section in manager.sections:
    print(f"Section with upper ID {section.upper} and lower ID {section.lower}")
