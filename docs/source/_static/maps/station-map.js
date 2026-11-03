// Leaflet maps for the site-lookup pages (bufkit_map, raob_map), adapted
// from the SounderPy web app's station-map.js. Needs the Leaflet and topojson-client
// globals; the pages load both before this file.

(() => {
  // Basemap drawn from Natural Earth data (public domain), so there are no tile servers or
  // API keys involved.
  const OCEAN = "#3a4749";
  const BASE = {
    land: { stroke: true, color: "#2c2820", weight: 1, fillColor: "#5c5443", fillOpacity: 1 },
    lakes: { stroke: true, color: "#2c2820", weight: 0.7, fillColor: OCEAN, fillOpacity: 1 },
    highways: { color: "#a86c4c", weight: 1, opacity: 0.75 },
    countries: { color: "#231f19", weight: 1.3, opacity: 0.9 },
    states: { color: "#38332a", weight: 0.8, opacity: 0.85 },
    graticule: { color: "#4a5556", weight: 1, dashArray: "2 4" },
  };
  const HIGHWAY_MIN_ZOOM = 5;
  const MAX_PLACE_LABELS = 120;
  const DATA = "_static/maps/data/";

  const DOT = { radius: 4.5, weight: 1.5, color: "#1a1712", fillColor: "#f3ead6", fillOpacity: 1 };
  const PICKED = { radius: 8, weight: 3, color: "#e0a93b", fillColor: "#ffffff", fillOpacity: 1 };

  // At most one station per CELL x CELL pixel square at the current zoom, so dense areas stay readable.
  const CELL = 44;

  const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => `&#${c.charCodeAt(0)};`);

  async function fetchData(file) {
    const res = await fetch(DATA + file);
    if (!res.ok) throw new Error(`${file}: HTTP ${res.status}`);
    return res.json();
  }

  async function loadBasemap() {
    const [world, us, lakes, highways, places] = await Promise.all(
      ["world-50m.json", "us-states-10m.json", "lakes-50m.json", "highways-na.json", "places.json"].map(fetchData),
    );
    return {
      land: topojson.feature(world, world.objects.land),
      countries: topojson.mesh(world, world.objects.countries, (a, b) => a !== b),
      states: topojson.mesh(us, us.objects.states, (a, b) => a !== b),
      lakes, highways, places,
    };
  }

  // Town labels: the most important places that fit the current zoom and view.
  function placeLabels(map, places) {
    const layer = L.layerGroup().addTo(map);
    const update = () => {
      const zoom = map.getZoom(), bounds = map.getBounds().pad(0.05);
      layer.clearLayers();
      let shown = 0;
      for (const [name, lat, lon, minZoom] of places) {
        if (minZoom > zoom || shown >= MAX_PLACE_LABELS) break;
        if (!bounds.contains([lat, lon])) continue;
        L.marker([lat, lon], {
          pane: "labels", interactive: false, keyboard: false,
          icon: L.divIcon({ className: "spy-place", html: `<span>${esc(name)}</span>`, iconSize: null }),
        }).addTo(layer);
        shown++;
      }
    };
    map.on("moveend zoomend", update);
    update();
  }

  function graticule(step = 10) {
    const lines = [];
    for (let lat = -80; lat <= 80; lat += step) lines.push([[lat, -180], [lat, 180]]);
    for (let lon = -180; lon <= 180; lon += step) lines.push([[-85, lon], [85, lon]]);
    return lines;
  }

  async function drawBasemap(map) {
    map.createPane("basemap").style.zIndex = 200;
    map.createPane("labels").style.zIndex = 350; // above the basemap, below station dots
    const renderer = L.canvas({ pane: "basemap" });
    const layer = { pane: "basemap", renderer, interactive: false };
    L.polyline(graticule(), { ...layer, ...BASE.graticule }).addTo(map);

    const { land, lakes, highways, countries, states, places } = await loadBasemap();
    for (const [shape, style] of [[land, BASE.land], [lakes, BASE.lakes]]) {
      L.geoJSON(shape, { ...layer, style }).addTo(map);
    }
    const roads = L.geoJSON(highways, { ...layer, style: BASE.highways }).addTo(map);
    for (const [shape, style] of [[countries, BASE.countries], [states, BASE.states]]) {
      L.geoJSON(shape, { ...layer, style }).addTo(map);
    }

    // Highways only once zoomed in to a region; at continent scale they're clutter.
    const toggleRoads = () => roads.setStyle({ opacity: map.getZoom() >= HIGHWAY_MIN_ZOOM ? BASE.highways.opacity : 0 });
    map.on("zoomend", toggleRoads);
    toggleRoads();
    placeLabels(map, places);
  }

  function createBaseMap(el) {
    const map = L.map(el, {
      center: [39, -96], zoom: 4, minZoom: 2, maxZoom: 9, preferCanvas: true,
      // On phones one finger scrolls the page; two fingers move the map.
      dragging: !L.Browser.mobile,
      zoomSnap: 0.5, wheelPxPerZoomLevel: 90,
      maxBounds: [[-85, -200], [85, 200]], maxBoundsViscosity: 0.8,
    });
    map.attributionControl.addAttribution('Map data: <a href="https://www.naturalearthdata.com/">Natural Earth</a>');
    drawBasemap(map).catch((err) => console.error("Basemap failed to load", err));
    if (L.Browser.mobile) {
      const note = L.control({ position: "bottomleft" });
      note.onAdd = () => Object.assign(L.DomUtil.create("div", "spy-map-touch-note"), { textContent: "Use two fingers to move the map" });
      note.addTo(map);
    }
    new ResizeObserver(() => map.invalidateSize()).observe(el);
    return map;
  }

  /**
   * Clickable station dots with faded ID labels, thinned out by zoom.
   * rows: [id, name, region, lat, lon, elev].
   */
  function stationLayer(map, rows) {
    const pane = map.createPane("station-labels");
    pane.style.zIndex = 450;
    pane.style.pointerEvents = "none";

    const entries = [];
    let picked = null;

    function pick(entry) {
      picked?.marker.setStyle(DOT);
      picked = entry;
      entry.marker.setStyle(PICKED);
      update();
    }

    // Main sites (4-character IDs like KOAX) win when stations crowd the same spot.
    const ordered = [...rows].sort((a, b) => (b[0].length === 4) - (a[0].length === 4));
    for (const row of ordered) {
      const [id, name, region, lat, lon, elev] = row;
      const entry = { row };
      entry.marker = L.circleMarker([lat, lon], DOT)
        .bindTooltip(`${esc(id)} · ${esc(name)}`, { direction: "top", offset: [0, -4] })
        .bindPopup(`<strong>${esc(id)}</strong><span>${esc([name, region].filter(Boolean).join(", "))}</span><small>${esc(elev)} m elevation</small>`, { className: "spy-popup" })
        .on("click", () => pick(entry));
      entry.label = L.marker([lat, lon], {
        pane: "station-labels", interactive: false, keyboard: false,
        icon: L.divIcon({ className: "spy-station-label", html: `<span>${esc(id)}</span>`, iconSize: null }),
      });
      entries.push(entry);
    }

    function update() {
      const zoom = map.getZoom(), bounds = map.getBounds().pad(0.2);
      const taken = new Set(), visible = new Set();
      const claim = (entry) => {
        const p = map.project(entry.marker.getLatLng(), zoom);
        const cell = `${Math.floor(p.x / CELL)}:${Math.floor(p.y / CELL)}`;
        if (taken.has(cell)) return;
        taken.add(cell);
        visible.add(entry);
      };
      if (picked) { claim(picked); visible.add(picked); }
      entries.filter((entry) => entry !== picked && bounds.contains(entry.marker.getLatLng())).forEach(claim);
      for (const entry of entries) {
        const show = visible.has(entry);
        for (const layer of [entry.marker, entry.label]) {
          if (show && !map.hasLayer(layer)) layer.addTo(map);
          else if (!show && map.hasLayer(layer)) layer.remove();
        }
      }
      picked?.marker.bringToFront();
    }
    map.on("moveend zoomend", update);
    update();
  }

  window.spyStationMap = async (source) => {
    const map = createBaseMap(document.getElementById("map"));
    stationLayer(map, await fetchData(`${source}-stations.json`));
  };
})();
