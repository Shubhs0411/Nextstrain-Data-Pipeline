"""
Location-to-country mapping for inferring country/region from strain names.
Used when BV-BRC metadata has empty Isolation Country.
"""
from typing import Optional, Tuple

# US states, territories, regions -> USA
US_LOCATIONS = {
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut",
    "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
    "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan",
    "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New_Hampshire",
    "New_Jersey", "New_Mexico", "New_York", "New_York_City", "North_Carolina", "North_Dakota",
    "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode_Island", "South_Carolina",
    "South_Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington",
    "WA", "West_Virginia", "Wisconsin", "Wyoming", "District_of_Columbia", "DC",
    "King_County", "NY", "LA", "SF",
}

# Cities/regions that map to country (location -> (country, region))
LOCATION_TO_COUNTRY: dict[str, Tuple[str, str]] = {
    # USA (country and states)
    "USA": ("USA", "North America"),
    "United_States": ("USA", "North America"),
    **{loc: ("USA", "North America") for loc in US_LOCATIONS},
    # China
    "Wuhan": ("China", "Asia"),
    "Nanjing": ("China", "Asia"),
    "Guangdong": ("China", "Asia"),
    "Shanghai": ("China", "Asia"),
    "Beijing": ("China", "Asia"),
    "Hong_Kong": ("China", "Asia"),
    "Hong Kong": ("China", "Asia"),
    "Zhejiang": ("China", "Asia"),
    "Hebei": ("China", "Asia"),
    "Foshan": ("China", "Asia"),
    "Gaoyou": ("China", "Asia"),
    "Ganzhou": ("China", "Asia"),
    "Liaoning": ("China", "Asia"),
    "Fujian": ("China", "Asia"),
    # Australia
    "Australia": ("Australia", "Oceania"),
    "Western_Australia": ("Australia", "Oceania"),
    "South_Australia": ("Australia", "Oceania"),
    "Victoria": ("Australia", "Oceania"),
    "Queensland": ("Australia", "Oceania"),
    "New_South_Wales": ("Australia", "Oceania"),
    "Brisbane": ("Australia", "Oceania"),
    "Sydney": ("Australia", "Oceania"),
    "Melbourne": ("Australia", "Oceania"),
    "Newcastle": ("Australia", "Oceania"),
    # New Zealand
    "New_Zealand": ("New Zealand", "Oceania"),
    "Christchurch": ("New Zealand", "Oceania"),
    "Auckland": ("New Zealand", "Oceania"),
    # Chile
    "Santiago": ("Chile", "South America"),
    "Chile": ("Chile", "South America"),
    # Other countries
    "England": ("United Kingdom", "Europe"),
    "Scotland": ("United Kingdom", "Europe"),
    "Wales": ("United Kingdom", "Europe"),
    "UK": ("United Kingdom", "Europe"),
    "Denmark": ("Denmark", "Europe"),
    "France": ("France", "Europe"),
    "Germany": ("Germany", "Europe"),
    "Italy": ("Italy", "Europe"),
    "Spain": ("Spain", "Europe"),
    "Netherlands": ("Netherlands", "Europe"),
    "Sweden": ("Sweden", "Europe"),
    "Norway": ("Norway", "Europe"),
    "Japan": ("Japan", "Asia"),
    "Korea": ("South Korea", "Asia"),
    "South_Korea": ("South Korea", "Asia"),
    "Thailand": ("Thailand", "Asia"),
    "Vietnam": ("Vietnam", "Asia"),
    "Singapore": ("Singapore", "Asia"),
    "India": ("India", "Asia"),
    "Bangladesh": ("Bangladesh", "Asia"),
    "Indonesia": ("Indonesia", "Asia"),
    "Malaysia": ("Malaysia", "Asia"),
    "Philippines": ("Philippines", "Asia"),
    "Cambodia": ("Cambodia", "Asia"),
    "Iran": ("Iran", "Asia"),
    "Israel": ("Israel", "Asia"),
    "Turkey": ("Turkey", "Asia"),
    "Russia": ("Russia", "Europe"),
    "Moscow": ("Russia", "Europe"),
    "Mongolia": ("Mongolia", "Asia"),
    "Egypt": ("Egypt", "Africa"),
    "Uganda": ("Uganda", "Africa"),
    "Kenya": ("Kenya", "Africa"),
    "South_Africa": ("South Africa", "Africa"),
    "Nigeria": ("Nigeria", "Africa"),
    "Canada": ("Canada", "North America"),
    "Manitoba": ("Canada", "North America"),
    "Ontario": ("Canada", "North America"),
    "Quebec": ("Canada", "North America"),
    "Mexico": ("Mexico", "North America"),
    "Monterrey": ("Mexico", "North America"),
    "Hidalgo": ("Mexico", "North America"),
    "Brazil": ("Brazil", "South America"),
    "Argentina": ("Argentina", "South America"),
    "Colombia": ("Colombia", "South America"),
    "Peru": ("Peru", "South America"),
    "Ecuador": ("Ecuador", "South America"),
    "Venezuela": ("Venezuela", "South America"),
    "Nicaragua": ("Nicaragua", "North America"),
    "Honduras": ("Honduras", "North America"),
    "El_Salvador": ("El Salvador", "North America"),
    "Guatemala": ("Guatemala", "North America"),
    "Panama": ("Panama", "North America"),
    "Costa_Rica": ("Costa Rica", "North America"),
    "Cuba": ("Cuba", "North America"),
    "Puerto_Rico": ("Puerto Rico", "North America"),
    "Jamaica": ("Jamaica", "North America"),
    "Dominican_Republic": ("Dominican Republic", "North America"),
    "Haiti": ("Haiti", "North America"),
    "Greece": ("Greece", "Europe"),
    "GREECE": ("Greece", "Europe"),
    "Dhaka": ("Bangladesh", "Asia"),
    "Sukhbaatar": ("Mongolia", "Asia"),
    "Ansan": ("South Korea", "Asia"),
    "Fukui": ("Japan", "Asia"),
    "Hokkaido": ("Japan", "Asia"),
    "Bethesda": ("USA", "North America"),
    "Memphis": ("USA", "North America"),
    "Gainesville": ("USA", "North America"),
    "Albany": ("USA", "North America"),
    "Alchi": ("Japan", "Asia"),
    "Milano": ("Italy", "Europe"),
    "Saraburi": ("Thailand", "Asia"),
    "CHL": ("Chile", "South America"),
}

# Country -> region for filling when region is empty
COUNTRY_TO_REGION: dict[str, str] = {
    "USA": "North America",
    "United States": "North America",
    "Canada": "North America",
    "Mexico": "North America",
    "China": "Asia",
    "Japan": "Asia",
    "South Korea": "Asia",
    "Thailand": "Asia",
    "Vietnam": "Asia",
    "Singapore": "Asia",
    "India": "Asia",
    "Bangladesh": "Asia",
    "Indonesia": "Asia",
    "Malaysia": "Asia",
    "Philippines": "Asia",
    "Cambodia": "Asia",
    "Iran": "Asia",
    "Australia": "Oceania",
    "New Zealand": "Oceania",
    "United Kingdom": "Europe",
    "UK": "Europe",
    "Denmark": "Europe",
    "France": "Europe",
    "Germany": "Europe",
    "Italy": "Europe",
    "Spain": "Europe",
    "Netherlands": "Europe",
    "Sweden": "Europe",
    "Norway": "Europe",
    "Greece": "Europe",
    "Russia": "Europe",
    "Mongolia": "Asia",
    "Turkey": "Asia",
    "Egypt": "Africa",
    "Uganda": "Africa",
    "Kenya": "Africa",
    "South Africa": "Africa",
    "Nigeria": "Africa",
    "Brazil": "South America",
    "Argentina": "South America",
    "Chile": "South America",
    "Colombia": "South America",
    "Peru": "South America",
    "Ecuador": "South America",
    "Venezuela": "South America",
    "Nicaragua": "North America",
    "Panama": "North America",
    "Cuba": "North America",
    "Puerto Rico": "North America",
    "Dominican Republic": "North America",
    "Haiti": "North America",
}


def infer_country_region_from_strain(strain: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse strain name (e.g. A/Human/Connecticut/Flu119/2013) and infer country/region.
    Returns (country, region) or (None, None) if not found.
    """
    parts = strain.replace(" ", "_").strip().split("/")
    if len(parts) < 3:
        return None, None
    # Skip A/ and host (e.g. Human)
    for i in range(2, len(parts) - 1):  # location is typically 3rd or 4th part
        loc = parts[i].strip()
        if loc and loc in LOCATION_TO_COUNTRY:
            country, region = LOCATION_TO_COUNTRY[loc]
            return country, region
    return None, None


def get_region_for_country(country: str) -> Optional[str]:
    """Return region for a given country when region is empty."""
    return COUNTRY_TO_REGION.get(country.strip(), None)
